#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$Bin/..";
use Parsing_Routines;
use Bad_SNPS_and_CNVS;
use Cwd 'realpath';
use Cwd;
use Getopt::Long;
use strict;

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $parsing = TCGA_Lib::Parsing_Routines->new;
my $Bad_SNP_CNV = TCGA_Lib::Bad_SNPS_and_CNVS->new;
my $Analysispath = realpath("../../Analysis");

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'output|o=s' => \my $bad_snps,#directory where the bad snps will be located
    'snp_dir|s=s' =>\my $snp_dir,#directories where the birdseed and copy number variants were downloaded(separate them by a comma)
                                 #The directory with the birdseed data downloaded mus be entered first otherwise this script will not run
    'affy_dir|a=s' => \my $affy_dir,
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("2.2");

my @raw_file_dirs;
if(defined $snp_dir)
{
    @raw_file_dirs = split(",",$snp_dir);
    if(scalar(@raw_file_dirs) < 2 or scalar(@raw_file_dirs) >= 3)
    {
        print "Enter in birdseed(Genotypes) and copy number variation(Copy number estimate) directories only.\n";
        $parsing->usage("2.2");
    }
}
else
{
    print "SNP directories were not entered\n";
    $parsing->usage("2.2");
}

if($help)
{
    $parsing->usage("2.2");
}

if(!defined $disease_abbr)
{
    print "disease type was not entered!\n";
    $parsing->usage("2.2");
}

my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";

if (!(-d $RNA_Path))
{
    print "$RNA_Path does not exist. Either it was deleted, moved or the scripts required for it have not been run.\n";
    exit;
}

if(!defined $bad_snps)
{
    $bad_snps = "$RNA_Path/bad_snps";
    print STDERR "using $bad_snps as default\n"; 
}

$affy_dir = "affy6" unless defined $affy_dir;
if(!defined $affy_dir)
{
    print "Enter in the the directory that was specified in script 1.0\n";
    $parsing->usage("2.2"); 
}

my $SNP_DIR1="$Analysispath/$disease_abbr/SNP6/$raw_file_dirs[0]";
my $SNP_DIR2="$Analysispath/$disease_abbr/SNP6/$raw_file_dirs[1]";

chdir "$RNA_Path";

`mkdir -p $bad_snps` unless(-d "$bad_snps");
`rm -f $RNA_Path/bad_snps/*`;

unless($bad_snps = "$RNA_Path/bad_snps")
{
   $bad_snps = realpath("../$bad_snps"); 
}

#make tumor normal tables for birdseed and copy number
#mk_tn_tables(directory where Genotype data was downloaded,raw file name(same as directory name but will be used for files))
$Bad_SNP_CNV->mk_tn_tables("$SNP_DIR1",$raw_file_dirs[0]);
#mk_tn_tables(directory where Copy number estimate data was downloaded,raw file name(same as directory name but will be used for files))
$Bad_SNP_CNV->mk_tn_tables("$SNP_DIR2",$raw_file_dirs[1]);

$parsing->vlookup("look.have_cd",3,"$raw_file_dirs[0]\_T",3,"1,2,3,4","y","$raw_file_dirs[0]\_1");
$parsing->vlookup("$raw_file_dirs[0]\_1",3,"$raw_file_dirs[0]\_N",3,2,"y","$raw_file_dirs[0]\_2");
$parsing->vlookup("$raw_file_dirs[0]\_2",3,"$raw_file_dirs[1]\_T",3,2,"y","$raw_file_dirs[1]\_1");
$parsing->vlookup("$raw_file_dirs[1]\_1",3,"$raw_file_dirs[1]\_N",3,2,"y","$raw_file_dirs[1]\_2");

$parsing->pull_column("$raw_file_dirs[1]\_2","5,6,9,3,8,10,1","lookup_TN");

`grep NaN -v lookup_TN |sort -k 3,3 > $disease_abbr\_lookup_TN`;

#run_snps will link genotype from tumor and normal samples
#It will call the suproutine flag_snps to find out bad snps
#These bad snps will be written into the files within the bad snp directory that was specified or the default one if a directory was not
print "Now running run_snps\n";
#run_snps(directory where Genotype data was downloaded,directory where Copy number estimate data was downloaded,TN file generated after the vlookups and pull_column,user defined directory from command line or default directory)
$Bad_SNP_CNV->run_snps("$SNP_DIR1,$SNP_DIR2","$RNA_Path/$disease_abbr\_lookup_TN","$bad_snps");

#get_cd_bed will try to lookup these bed snps with affy6.cd.txt
#The results will be saved into the directory bad_snps_bed
print "Now running get_cd_bed.\n";
mkdir "$RNA_Path/bad_snps_bed";
`rm -f $RNA_Path/bad_snps_bed/*`;
#get_cd_bed(path to snp6.cd.txt,path to the RNA_Seq_Analysis directory,user defined directory from command line or default directory)
$Bad_SNP_CNV->get_cd_bed("$RNA_Path/$affy_dir/snp6.cd.txt","$RNA_Path","$bad_snps");

#run_cnv will use vlookup, strip_head, smooth, and delin_cnv.
# As CNV only have signal values but not genotypes, this will get the mean difference
# of signals between whatever number of CNVs spreading over the middle position of the compared CNV.
# smooth called from run_cnv is defined for this job.
# Its output would be chr middle_pos mean_sig_of_#_CNVs

# run_cnv also calls delin_cnv, which is going to keep the start and end of a zone with
# all consecutive mean_sig_of_#_CNV > 0.5
# This can be easily virtualized in the following plots
# Only the two consecutive positions of '?' will be kept as
# the START and END positions for Bad CNVs!
#                |
#            0.5 |                     ?--------------?
# signal diff    |                     
#                |---------------------                ----------------------
#                |
#                ----------------------------------------------------------------|
#            
#                |
#            0.5 |?-----------?                    ?------------------------------------?
# signal diff    |
#                |             -------------------
#                |
#                ----------------------------------------------------------------|
#                          chromosome position

#bad_cnv.bed files will be saved into the dir bad_cnvs
#If there are errors in this process then it most likely means that some CNVs were not download properly.
print "Now running run_cnv.\n";
mkdir "$RNA_Path/bad_cnvs";
`rm -f $RNA_Path/bad_cnvs/*`;
#run_cnv(TN file generated after the vlookups and pull_column,directory where Copy number estimate data was downloaded,directory where Genotype data was downloaded,cnv.hg19.bed file,bad_cnvs directory)
$Bad_SNP_CNV->run_cnv("$RNA_Path/$disease_abbr\_lookup_TN","$SNP_DIR2","$SNP_DIR1","$RNA_Path/$affy_dir/cnv.hg19.bed","$RNA_Path/bad_cnvs");

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;