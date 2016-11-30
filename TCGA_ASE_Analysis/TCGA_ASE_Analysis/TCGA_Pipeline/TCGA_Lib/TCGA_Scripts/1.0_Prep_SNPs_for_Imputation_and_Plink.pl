#!/usr/bin/perl -w

use MCE;
use FindBin qw($Bin);
use lib "$Bin/..";
use Imputation_Plink;
use Parsing_Routines;
use Cwd 'realpath';
use Getopt::Long;
use strict;

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $imput_plink = TCGA_Lib::Imputation_Plink->new;
my $parsing = TCGA_Lib::Parsing_Routines->new;
my $TCGA_Pipeline_Dir = realpath("../../");

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'output|o=s' => \my $affy_dir,
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage;

if ($help)
{
    $parsing->usage;
}

if (!defined $disease_abbr)
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

my $database_path = "$TCGA_Pipeline_Dir/Database";

#Checks if there is no Database directory
if(!(-d "$database_path"))
{
    print STDERR "$database_path does not exist, it was either moved, deleted or has not been downloaded.\nPlease check the README.md file on the github page to find out where to get the Database directory.\n";
    exit;
}

my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";
mkdir "$Analysispath/$disease_abbr/RNA_Seq_Analysis" unless(-d "$Analysispath/$disease_abbr/RNA_Seq_Analysis");

if (!defined $affy_dir)
{
    $affy_dir = "affy6";
    print STDERR "using $affy_dir as default\n"; 
}

chdir $RNA_Path;

$affy_dir =~ s/\/$//;
`mkdir -p $RNA_Path/$affy_dir` unless(-d "$RNA_Path/$affy_dir");

if (-e "$RNA_Path/$affy_dir/snp6.anno.txt" and -e "$RNA_Path/$affy_dir/snp6.cd.txt") 
{
    print "the necessary files (snp6.anno.txt and snp6.cd.txt) have already been processed and reside in $RNA_Path/$affy_dir/\n";
    print "proceed to script 1.1\n";
    exit;
}
else
{
    unless (-e "$database_path/GenomeWideSNP_6.na35.annot.csv")
    {
        `unzip $database_path/GenomeWideSNP_6.na35.annot.csv.zip -d $RNA_Path/$affy_dir`;
    }
    
    #ParseAnnoAffxSNP6 changes the strand '-' into "+" and transforms two alleles in /ATCG/ into their complements /TAGC/
    print "Now running ParseAnnoAffxSNP6\n";
    #ParseAnnoAffxSNP6(GenomeWideSNP_6.na35.annot.csv file, output file)
    $imput_plink->ParseAnnoAffxSNP6("$RNA_Path/$affy_dir/GenomeWideSNP_6.na35.annot.csv","$RNA_Path/$affy_dir/snp6.anno.txt");
    
    #Get snp6.cd.txt for imputation, which is saved in the dir affy6
    print "Now running PrepareAffxAlleleAccording2ImputationAlleles\n";
    #PrepareAffxAlleleAccording2ImputationAlleles(path to GenomeWideSNP_6.na35.annot.csv,database directory path)
    $imput_plink->PrepareAffxAlleleAccording2ImputationAlleles("$RNA_Path/$affy_dir","$database_path");
    
    print "All jobs have finished for $disease_abbr\n";
    
    $time = localtime;
    print "Script finished on $time.\n";
    
    exit;
}