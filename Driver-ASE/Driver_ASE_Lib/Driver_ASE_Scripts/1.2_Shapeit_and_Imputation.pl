#!/usr/bin/perl -w

use MCE::Map;
use FindBin qw($Bin);
use lib "$Bin/..";
use Imputation_Plink;
use Parsing_Routines;
use Cwd 'realpath';
use Cwd;
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
    'shapeit|s=s' => \my $shapeit,#path to shapeit
    'impute|i=s' => \my $impute,#path to impute2
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("1.2");

if($help)
{
    $parsing->usage("1.2");
}

if(!defined $cancer_type)
{
    print STDERR "Cancer type was not entered!\n";
    $parsing->usage("1.2");
}

if (!defined $shapeit)
{
    $shapeit = "shapeit";
}
if (!defined $impute)
{
    $impute = "impute2";
}

$parsing->CheckSoftware("$shapeit","$impute");

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $RNA_Path = "$Analysispath/$cancer_type/RNA_Seq_Analysis";
my $map_dir = "maps";
my $ped_dir = "peds";
my $logs = "logs";
my $imputation = "$RNA_Path/phased";
my $OneKG_Ref_Path = "$database_path/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono";
my $Impute2out = "$Analysispath/$cancer_type/phased_imputed_raw_out";

$parsing->check_directory_existence("$database_path","$OneKG_Ref_Path","$Analysispath","$Analysispath/$cancer_type","$RNA_Path","$RNA_Path/$ped_dir","$RNA_Path/$map_dir"); #check if directories or files exist
$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

chdir "$RNA_Path";

`mkdir -p $imputation` unless(-d "$imputation");
`rm -f $imputation/*`;
`mkdir -p $logs` unless(-d "$logs");

#submit_shapeit($OneKG_Ref_Path (path to ALL.integrated_phase1_SHAPEIT_16-06-14.nomono),$RNA_Path (path to directory where RNA-Seq analysis data is stored),$imputation (path to directory with default name phased),$map_dir (path to directory where map files are stored),$ped_dir(path to directory where ped files are stored),l$logs (path to directory where log files will be stored),$shapeit (path/command for shapeit))
my @shapeit_cmds = $impute_plink->submit_shapeit("$OneKG_Ref_Path","$RNA_Path","$imputation","$map_dir","$ped_dir","$logs","$shapeit");

mce_map
{
    system("$shapeit_cmds[$_]");
} 0..$#shapeit_cmds;

#fetch_Chrom_Sizes(reference genome(e.g. hg19))
$impute_plink->fetch_Chrom_Sizes("hg19");

`cat chr_lens_grep_chr | grep -v Un | grep -v random | grep -v hap | grep -v M | grep -v Y | grep chr > chr_lens`;

`cat chr_lens|head -n 11 > file_for_submit`;

mkdir "$Impute2out" unless(-d "$Impute2out");
#submit first 11 chrs for imputation
#submit_all(file with chr sizes,$OneKG_Ref_Path (path to ALL.integrated_phase1_SHAPEIT_16-06-14.nomono),$imputation (path to directory with default name phased),$Impute2out (path to directory where raw imputed files will be stored))
my @imput2cmds = $impute_plink->submit_all("file_for_submit", $OneKG_Ref_Path, $imputation, $Impute2out,"$impute");

mce_map
{
   system("$imput2cmds[$_]");
} 0..$#imput2cmds;

#To save disk space;
#remove them permentately!
`rm -f $Impute2out/*_allele_probs`;
`rm -f $Impute2out/*_info_by_sample`;
`rm -f $Impute2out/*_info`;
`rm -f $Impute2out/*_summary`;
`rm -f $Impute2out/*_warnings`;

`cat chr_lens|tail -n 12 > file_for_submit`;
#submit next 12 chrs for imputation
undef @imput2cmds;
@imput2cmds = $impute_plink->submit_all( "file_for_submit", $OneKG_Ref_Path, $imputation, $Impute2out,"$impute");
mce_map
{
    system("$imput2cmds[$_]");
} 0..$#imput2cmds;

#To save disk space;
#remove them permentately!
`rm -f $Impute2out/*_allele_probs`;
`rm -f $Impute2out/*_info_by_sample`;
`rm -f $Impute2out/*_info`;
`rm -f $Impute2out/*_summary`;
`rm -f $Impute2out/*_warnings`;

print "All jobs are done for $cancer_type.\n";

$time = localtime;
print "Script finished: $time.\n";

exit;
