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
no warnings 'once';

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

if($help)
{
    $parsing->usage;
}

if(!defined $cancer_type)
{
    print STDERR "Cancer type was not entered!\n";
    $parsing->usage;
}

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $RNA_Path = "$Analysispath/$cancer_type/RNA_Seq_Analysis";
my $affy_dir = "affy6";
my $SNP6_Raw_Files_Dir = "$Analysispath/$cancer_type/SNP6/Copy number estimate";
my $GenomeWideSNP = "GenomeWideSNP_6.cn.na35.annot.csv.zip";
my $cds_sorted = "cds_sorted";
my $cnv_hg19_bed = "cnv.hg19.bed";
my $look_have_cd = "look.have_cd";

$parsing->check_directory_existence("$database_path","$database_path/$GenomeWideSNP","$Analysispath","$Analysispath/$cancer_type","$RNA_Path","$RNA_Path/$affy_dir","$RNA_Path/$cds_sorted","$SNP6_Raw_Files_Dir"); 
#check if directories or files exist
$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

chdir "$RNA_Path";

print "Running zcat on $GenomeWideSNP.\n";
open (F,"zcat $database_path/$GenomeWideSNP|");
open (O,">annot_csv");
while(my $r = <F>)
{
  print O $r;      
}
close (F);
close (O);

print "Running pull_bed\n";
#pull_bed($RNA_Path/annot_csv (annot_csv file create from $GenomeWideSNP in the $database_path directory), $RNA_Path/$affy_dir/$cnv_hg19_bed (path to cnv.hg19.bed))
$Bad_SNP_CNV->pull_bed("$RNA_Path/annot_csv","$RNA_Path/$affy_dir/$cnv_hg19_bed");

`ls "$SNP6_Raw_Files_Dir" > cnv_file_names.txt`;

#pull normal CNV samples
print "Now running pull_normal.\n";
#pull_normal(file containing a list of cnv files,output file)
$Bad_SNP_CNV->pull_normal("cnv_file_names.txt","normal_cnv.txt");

#move barcode from the cds_sorted directory into a file called have_cd
print "Now running mv_bcode for cds_sorted_dir.\n";
`ls $RNA_Path/$cds_sorted > cds_sorted_dir`;
#mv_bcode(file containing list of bed files from the directory with sorted cds files,output file)
$Bad_SNP_CNV->mv_bcode("cds_sorted_dir","have_cd");

print "Now running hv_cd.\n";
$parsing->vlookup("normal_cnv.txt",2,"have_cd",2,1,"y","have_hv_cd");
#hv_cd(output file from vlookup,output file)
$Bad_SNP_CNV->hv_cd("have_hv_cd","hv_cd_out");
print "Now running dump_non_01_10_11.\n";
#dump_non_01_10_11(output file from hv_cd,output file)
$Bad_SNP_CNV->dump_non_01_10_11("hv_cd_out","$look_have_cd");

print "All jobs have finished for $cancer_type.\n";

$time = localtime;
print "Script finished: $time.\n";

exit;
