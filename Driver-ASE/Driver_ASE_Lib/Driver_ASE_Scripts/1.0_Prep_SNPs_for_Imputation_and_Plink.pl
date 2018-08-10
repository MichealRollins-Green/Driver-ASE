#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$Bin/..";
use Imputation_Plink;
use Parsing_Routines;
use Cwd 'realpath';
use Getopt::Long;
use strict;
use autodie;

my $time = localtime;
print "Script started: $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $impute_plink = Driver_ASE_Lib::Imputation_Plink->new;
my $parsing = Driver_ASE_Lib::Parsing_Routines->new;

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
my $Analysispath = realpath("../../Analysis");
my $database_path = "$Driver_ASE_Dir/Database";
my $RNA_Path = "$Analysispath/$cancer_type/RNA_Seq_Analysis";
my $affy_dir = "affy6";
my $snp6_anno_file = "snp6.anno.txt";
my $snp6_cd_file = "snp6.cd.txt";
my $GenomeWideSNP = "GenomeWideSNP_6.na35.annot.csv";

$parsing->check_directory_existence("$database_path","$Analysispath/$cancer_type"); #check if directories or files exist
$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

mkdir "$RNA_Path" unless(-d "$RNA_Path");

chdir $RNA_Path;

$affy_dir =~ s/\/$//;
`mkdir -p $RNA_Path/$affy_dir` unless(-d "$RNA_Path/$affy_dir");

if ((-e "$RNA_Path/$affy_dir/$snp6_anno_file" and -s "$RNA_Path/$affy_dir/$snp6_anno_file" !=0 ) and (-e "$RNA_Path/$affy_dir/$snp6_cd_file" and -s "$RNA_Path/$affy_dir/$snp6_cd_file" != 0)) 
{
    print "The necessary files ($snp6_anno_file and $snp6_cd_file) have already been processed and reside in $RNA_Path/$affy_dir/\n";
    print "proceed to script 1.1_Birdseed_to_ped_and_maps.pl\n";
}
else
{
    unless (-e "$RNA_Path/$affy_dir/$GenomeWideSNP")
    {
        print "unzipping $GenomeWideSNP\n";
        `unzip $database_path/$GenomeWideSNP.zip -d $RNA_Path/$affy_dir`;
    }
    
    #ParseAnnoAffxSNP6 changes the strand '-' into "+" and transforms two alleles in /ATCG/ into their complements /TAGC/
    print "Now running ParseAnnoAffxSNP6.\n";
    #ParseAnnoAffxSNP6($RNA_Path/$affy_dir/$GenomeWideSNP (path to GenomeWideSNP_6.na35.annot.csv file), $RNA_Path/$affy_dir/$snp6_anno_file(output file))
    $impute_plink->ParseAnnoAffxSNP6("$RNA_Path/$affy_dir/$GenomeWideSNP","$RNA_Path/$affy_dir/$snp6_anno_file");
    
    #Get $snp6_cd_file for imputation, which is saved in the dir affy6
    print "Now running PrepareAffxAlleleAccording2ImputationAlleles.\n";
    #PrepareAffxAlleleAccording2ImputationAlleles($RNA_Path/$affy_dir/$GenomeWideSNP (path to GenomeWideSNP_6.na35.annot.csv file),$database_path (path to the Database directory), $snp6_anno_file (snp6.anno.txt file),$snp6_cd_file (snp6.cd.txt file))
    $impute_plink->PrepareAffxAlleleAccording2ImputationAlleles("$RNA_Path/$affy_dir","$database_path","$snp6_anno_file","$snp6_cd_file");
    
    print "All jobs have finished for $cancer_type.\n";
    
    $time = localtime;
    print "Script finished: $time.\n";
}

exit;
