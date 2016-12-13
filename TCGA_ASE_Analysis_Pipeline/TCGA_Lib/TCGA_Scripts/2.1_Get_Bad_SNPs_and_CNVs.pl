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

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage;

if($help)
{
    $parsing->usage;
}

if(!defined $disease_abbr)
{
    print "disease type was not entered!\n";
    $parsing->usage;
}

my $Analysispath = realpath("../../Analysis");

#Checks if there is no Analysis directory
if (!(-d "$Analysispath"))
{
    print STDERR "$Analysispath does not exist, it was either deleted, moved or the script that creates it wasn't ran.\n";
    exit;
}
elsif(!(-d "$Analysispath/$disease_abbr"))
{
    print STDERR "$Analysispath/$disease_abbr does not exist, it was either deleted, moved or the script that creates it wasn't ran.\n";
    exit;
}

my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";

if (!(-d $RNA_Path))
{
    print "$RNA_Path does not exist. Either it was deleted, moved or the scripts required for it have not been run.\n";
    exit;
}

my $bad_snps = "$RNA_Path/bad_snps";

my $affy_dir = "affy6";
if(!(-d "$RNA_Path/$affy_dir"))
{
    print STDERR "The directory $RNA_Path/$affy_dir does not exist, it was moved, renamed. deleted or the script that creates it has not run.\n";
    exit; 
}

my $geno_dir = "$Analysispath/$disease_abbr/SNP6/Genotypes";
my $copy_dir = "$Analysispath/$disease_abbr/SNP6/Copy number estimate";

chdir "$RNA_Path";

`mkdir -p $bad_snps` unless(-d "$bad_snps");
`rm -f $RNA_Path/$bad_snps/*`;

#make tumor normal tables for birdseed and copy number
#mk_tn_tables(directory where Genotype data was downloaded,raw file name(same as directory name but will be used for files))
$Bad_SNP_CNV->mk_tn_tables("$geno_dir","Genotypes");
#mk_tn_tables(directory where Copy number estimate data was downloaded,raw file name(same as directory name but will be used for files))
$Bad_SNP_CNV->mk_tn_tables("$copy_dir","Copy number estimate");

$parsing->vlookup("look.have_cd",3,"Genotypes_T",3,"1,2,3,4","y","Genotypes_1");
$parsing->vlookup("Genotypes_1",3,"Genotypes_N",3,2,"y","Genotypes_2");
$parsing->vlookup("Genotypes_2",3,"Copy number estimate_T",3,2,"y","Copy number estimate_1");
$parsing->vlookup("Copy number estimate_1",3,"Copy number estimate_N",3,2,"y","Copy number estimate_2");

$parsing->pull_column("Copy number estimate_2","5,6,9,3,8,10,1","lookup_TN");

`grep NaN -v lookup_TN |sort -k 3,3 > $disease_abbr\_lookup_TN`;

#run_snps will link genotype from tumor and normal samples
#It will call the suproutine flag_snps to find out bad snps
#These bad snps will be written into the files within the bad snp directory that was specified or the default one if a directory was not
print "Now running run_snps\n";
#run_snps(directory where Genotype data was downloaded,directory where Copy number estimate data was downloaded,TN file generated after the vlookups and pull_column,user defined directory from command line or default directory)
$Bad_SNP_CNV->run_snps("$geno_dir,$copy_dir","$RNA_Path/$disease_abbr\_lookup_TN","$bad_snps");

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
$Bad_SNP_CNV->run_cnv("$RNA_Path/$disease_abbr\_lookup_TN","$copy_dir","$geno_dir","$RNA_Path/$affy_dir/cnv.hg19.bed","$RNA_Path/bad_cnvs");

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;