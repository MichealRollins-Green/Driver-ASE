#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$Bin/..";
use Parsing_Routines;
use Bad_SNPS_and_CNVS;
use Cwd 'realpath';
use Cwd;
use Getopt::Long;
use strict;
use autodie;

my $time = localtime;
print "Script started: $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $parsing = Driver_ASE_Lib::Parsing_Routines->new;
my $Bad_SNP_CNV = Driver_ASE_Lib::Bad_SNPS_and_CNVS->new;

GetOptions(
    'cancer|c=s' => \my $cancer_type,#e.g. OV
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage;

if ($help)
{
    $parsing->usage;
}

if (!defined $cancer_type)
{
    print STDERR "Cancer type was not entered!\n";
    $parsing->usage;
}

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $RNA_Path = "$Analysispath/$cancer_type/RNA_Seq_Analysis";
my $bad_snps = "$RNA_Path/bad_snps";
my $affy_dir = "affy6";
my $geno_dir = "$Analysispath/$cancer_type/SNP6/Genotypes";
my $copy_dir = "$Analysispath/$cancer_type/SNP6/Copy number estimate";
my $genotype = "Genotypes";
my $cnv = "Copy number estimate";
my $bad_snps_bed = "$RNA_Path/bad_snps_bed";
my $bad_cnvs = "bad_cnvs";
my $temp_dir = "temp";
my $snp6_cd_file = "snp6.cd.txt";
my $cnv_hg19_bed = "cnv.hg19.bed";
my $look_have_cd = "look.have_cd";

$parsing->check_directory_existence("$database_path","$Analysispath","$Analysispath/$cancer_type","$RNA_Path","$RNA_Path/$affy_dir","$geno_dir","$copy_dir"); #check if directories or files exist
$parsing->check_cancer_type($database_path,$cancer_type);

chdir "$RNA_Path";

`mkdir -p $bad_snps` unless(-d "$bad_snps");
`find $bad_snps/* 2>/dev/null|xargs rm -rf`;
#`rm -rf $bad_snps/*`;

#make tumor normal tables for birdseed and copy number
#mk_tn_tables($geno_dir (directory where Genotype data was downloaded), $genotype (raw file name (same as directory name but will be used for files)))
$Bad_SNP_CNV->mk_tn_tables("$geno_dir","$genotype");
#mk_tn_tables($copy_dir (directory where Copy number estimate data was downloaded), $cnv (raw file name (same as directory name but will be used for files)))
$Bad_SNP_CNV->mk_tn_tables("$copy_dir","$cnv");

$parsing->vlookup("$look_have_cd",2,"$genotype\_T",3,"1,2,3,4","y","$genotype\_1");
$parsing->vlookup("$genotype\_1",2,"$genotype\_N",3,2,"y","$genotype\_2");
$parsing->vlookup("$genotype\_2",2,"$cnv\_T",3,2,"y","$cnv\_1");
$parsing->vlookup("$cnv\_1",2,"$cnv\_N",3,2,"y","$cnv\_2");

$parsing->pull_column("$cnv\_2","4,6,9,3,8,10,1","lookup_TN");

`grep NaN -v lookup_TN |sort -k 3,3 > $cancer_type\_lookup_TN`;

#run_snps will link genotype from tumor and normal samples
#It will call the suproutine flag_snps to find out bad snps
#These bad snps will be written into the files within the bad snp directory that was specified or the default one if a directory was not
print "Now running run_snps.\n";
#run_snps($geno_dir,$copy_dir (directory where Genotype data was downloaded, and directory where Copy number estimate data was downloaded),$RNA_Path/$cancer_type\_lookup_TN (TN file generated after the vlookups and pull_column), $bad_snps (directory with default name bad_snps))
$Bad_SNP_CNV->run_snps("$geno_dir,$copy_dir","$RNA_Path/$cancer_type\_lookup_TN","$bad_snps");

#get_cd_bed will try to lookup these bed snps with affy6.cd.txt
#The results will be saved into the directory bad_snps_bed
print "Now running get_cd_bed.\n";
mkdir "$bad_snps_bed" unless(-d "$bad_snps_bed");
`find $bad_snps_bed/* 2>/dev/null|xargs rm -rf`;
#`rm -rf $bad_snps_bed/*`;
#get_cd_bed($RNA_Path/$affy_dir/$snp6_cd_file (path to snp6.cd.txt,path to the RNA_Seq_Analysis directory), $bad_snps_bed (directory with default name bad_snps_bed), $bad_snps (directory with default name bad_snps))
$Bad_SNP_CNV->get_cd_bed("$RNA_Path/$affy_dir/$snp6_cd_file","$bad_snps_bed","$bad_snps");

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
mkdir "$bad_cnvs" unless(-d "$bad_cnvs");
`find $bad_cnvs/* 2>/dev/null|xargs rm -rf`;
#`rm -rf $bad_cnvs/*`;
`mkdir $temp_dir` unless(-d "$temp_dir");
`find $temp_dir/* 2>/dev/null|xargs rm -rf`;
#`rm -rf $temp_dir/*`;
#run_cnv($RNA_Path/$cancer_type\_lookup_TN (TN file generated after the vlookups and pull_column), $copy_dir","$geno_dir (directory where Copy number estimate data was downloaded, and directory where Genotype data was downloaded), $RNA_Path/$affy_dir/$cnv_hg19_bed (cnv.hg19.bed file), $RNA_Path/$bad_cnvs (directory with default name bad_cnvs), $temp_dir (path to a temp directory))
$Bad_SNP_CNV->run_cnv("$RNA_Path/$cancer_type\_lookup_TN","$copy_dir","$geno_dir","$RNA_Path/$affy_dir/$cnv_hg19_bed","$RNA_Path/$bad_cnvs","$temp_dir");

print "All jobs have finished for $cancer_type.\n";

$time = localtime;
print "Script finished: $time.\n";

exit;
