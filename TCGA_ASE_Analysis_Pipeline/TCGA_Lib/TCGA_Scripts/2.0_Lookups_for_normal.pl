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
my $TCGA_Pipeline_Dir = realpath("../../");

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

my $database_path = "$TCGA_Pipeline_Dir/Database";

#Checks if there is no Database directory
if(!(-d "$database_path"))
{
    print STDERR "$database_path does not exist. It was either moved, renamed, deleted or has not been downloaded.\nPlease check the README.md file on the github page to find out where to get the Database directory.\n";
    exit;
}

my $Analysispath = realpath("../../Analysis");

#Checks if there is no Analysis directory
if (!(-d "$Analysispath"))
{
    print STDERR "$Analysispath does not exist. It was either deleted, moved or renamed.\n";
    print STDERR "Please run script 0_Download_SNPArray_From_GDC.pl.\n";
    exit;
}
elsif(!(-d "$Analysispath/$disease_abbr"))
{
    print STDERR "The Directory $Analysispath/$disease_abbr does not exist. It was either deleted, moved or renamed.\n";
    print STDERR "Please run script 0_Download_SNPArray_From_GDC.pl.\n";
    exit;
}

my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";

if (!(-d $RNA_Path))
{
    print STDERR "$RNA_Path does not exist. Either it was deleted, moved or renamed.\n";
    print STDERR "Please run script 1.0_Prep_SNPs_for_Imputation_and_Plink.pl.\n";
    exit;
}

my $affy_dir = "affy6";
if(!(-d "$RNA_Path/$affy_dir"))
{
    print "The directory $RNA_Path/$affy_dir does not exist. It was moved, renamed or deleted.\n";
    print "Please run script 1.0_Prep_SNPs_for_Imputation_and_Plink.pl.\n";
    exit; 
}

my $SNP6_Raw_Files_Dir = "$Analysispath/$disease_abbr/SNP6/Copy number estimate";

if(!(-d "$SNP6_Raw_Files_Dir"))
{
    print STDERR "The directory $SNP6_Raw_Files_Dir does not exist. It was moved, renamed or deleted";
    print STDERR "Please run script 0_Download_SNPArray_From_GDC.pl.\n";
    exit;
}

chdir "$RNA_Path";

print "Running zcat on GenomeWideSNP_6.cn.na35.annot.csv.zip\n";
open(F,"zcat $database_path/GenomeWideSNP_6.cn.na35.annot.csv.zip|") or die("Can't open zip for input: $!\n");
open(O,">annot_csv") or die("Can't open file for output");
while(my $r = <F>)
{
  print O $r;      
}
close(F);
close(O);

print "Running pull_bed\n";
#pull_bed(annot_csv file create from GenomeWideSNP_6.cn.na35.annot.csv.zip in the Database directory,path to cnv.hg19.bed)
$Bad_SNP_CNV->pull_bed("$RNA_Path/annot_csv","$RNA_Path/$affy_dir/cnv.hg19.bed");

`ls "$SNP6_Raw_Files_Dir" > cnv_file_names.txt`;

#pull normal CNV samples
print "Now running pull_normal\n";
#pull_normal(file containing a list of cnv files,output file)
$Bad_SNP_CNV->pull_normal("cnv_file_names.txt","normal_cnv.txt");

#move barcode from the cds_sorted directory into a file called have_cd
print "Now running mv_bcode for cds_sorted_dir\n";
`ls $RNA_Path/cds_sorted > cds_sorted_dir`;
#mv_bcode(file containing list of bed files from the cds_sorted directory,output file)
$Bad_SNP_CNV->mv_bcode("cds_sorted_dir","have_cd");

print "Now running hv_cd\n";
$parsing->vlookup("normal_cnv.txt",3,"have_cd",2,1,"y","have_hv_cd");
#hv_cd(output file from vlookup,output file)
$Bad_SNP_CNV->hv_cd("have_hv_cd","hv_cd_out");
print("Now running dump_non_01_10_11\n");
#dump_non_01_10_11(output file from hv_cd,output file)
$Bad_SNP_CNV->dump_non_01_10_11("hv_cd_out","look.have_cd");

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;