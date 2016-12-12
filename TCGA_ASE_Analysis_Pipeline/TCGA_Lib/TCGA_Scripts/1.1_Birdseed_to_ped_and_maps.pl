#!/usr/bin/perl -w

use MCE::Map;
use FindBin qw($Bin);
use lib "$Bin/..";
use Parsing_Routines;
use Imputation_Plink;
use Cwd 'realpath';
use Getopt::Long;
use strict;

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $parsing = TCGA_Lib::Parsing_Routines->new;
my $impute_plink = TCGA_Lib::Imputation_Plink->new;

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'output|o=s' => \my $bases_dir,
    'snp_dir|s=s' => \my $snp_dir,
    'affy_dir|a=s' => \my $affy_dir,
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("1.1");

if($help)
{
    $parsing->usage("1.1");
}

if(!defined $disease_abbr || !defined $snp_dir)
{
    print "disease type and/or snp directory was not entered!\n";
    $parsing->usage("1.1");
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

$affy_dir = "affy6" unless defined $affy_dir;
if(!(-d "$RNA_Path/$affy_dir"))
{
    print "Enter in the the directory that was specified in script 1.0 or run script 1.0 if it hasn't been run yet.\n";
    $parsing->usage("1.1"); 
}

chdir "$RNA_Path";

if(!defined $bases_dir)
{
    $bases_dir = "$RNA_Path/bases";
    print STDERR "using $bases_dir as default\n"; 
}

`mkdir -p $bases_dir` unless(-d "$bases_dir");
`rm -f $bases_dir/*`;

my $SNP6_Raw_Files_Dir="$Analysispath/$disease_abbr/SNP6/$snp_dir";

#gets all of the files in the the directory that was specified to download the SNP-Array data
my @birdseed = $parsing->get_only_files_in_dir("$SNP6_Raw_Files_Dir");

mce_map
{
    chomp ($_);
    #Prepare_Genos(Genotype files,directory where Genotype files were downloaded,snp6.cd.txt file created from 1.0,user or default output directory from command line,sample type(e.g. 0(normal), 1(tumor) or 2(normal/tumor))
    $impute_plink->Prepare_Genos("$_","$SNP6_Raw_Files_Dir","$RNA_Path/$affy_dir/snp6.cd.txt","$bases_dir",2);
}@birdseed;

`rm $RNA_Path/bases/*.t`;

my @files = $parsing->get_only_files_in_dir("$bases_dir");
   @files=grep{!/\.$/}@files;
my $selected_file=$files[-1];
#pull_column(file to pull column from,column(s) to pull, output file)
$parsing->pull_column("$bases_dir/$selected_file","1","$disease_abbr\_t");

mkdir "maps" unless(-d "$RNA_Path/maps");
`rm -f $RNA_Path/maps/*`;

$parsing->vlookup("$disease_abbr\_t",1,"$RNA_Path/$affy_dir/snp6.cd.txt",4,"1,3","y","map_for_sort");

`sort -k 2,2 -k 3,3n map_for_sort > map_for_pull_column`;
$parsing->pull_column("map_for_pull_column","2,3,1","mk_map_file");
#mk_map(file used to make maps, path to RNA_Seq_Analysis directory)
$impute_plink->mk_map("mk_map_file","$RNA_Path");

my @maps = $parsing->get_only_files_in_dir("$RNA_Path/maps");
@maps = grep{!/\.$/ and /^[0-9x]+.map$/i}@maps;
open(CHR,">chrs") or die("Can't open filehandle: $!");
for(my $i=0;$i<scalar(@maps);$i++)
{
    print CHR $maps[$i], "\n";
}
close(CHR);

mkdir "$RNA_Path/peds" unless(-d "$RNA_Path/peds");
`rm -f $RNA_Path/peds/*`;

#mk_ped_list_4_plink(file with chr.map,user or default output directory from command line,outfile for peds,path to RNA_Seq_Analysis directory)
$impute_plink->mk_ped_list_4_plink("chrs","$bases_dir","ped_list","$RNA_Path");

mce_map_f
{
    chomp($_);
    my @lines = split("\t",$_);
    #only keep normal samples for plink
    my $tn = [split("-",$lines[1])]->[-1];
    $tn =~ s/\D+//;
    #Make_Plink_Bed(chr#,sample(e.g. 10),sample id,path to RNA_Seq_Analysis directory,user or default output directory from command line)
    $impute_plink->Make_Plink_Bed("$lines[0]","$lines[1]","$lines[2]","$lines[3]","$bases_dir") if $tn > 9;
}"$RNA_Path/ped_list";

my @maps_bim = $parsing->get_only_files_in_dir("$RNA_Path/maps");
@maps_bim = grep {/bim/}@maps_bim;
open(BIMS,">$RNA_Path/maps/Small_Bims.txt") or die("Can't open file for output: $!");
for(my $i=0;$i<scalar(@maps_bim);$i++)
{
    $maps_bim[$i] =~ s/.bim//;
    print BIMS "$maps_bim[$i]\n";
}
close(BIMS);

chdir "$RNA_Path/maps";
#Merge all bims in dir maps with plink
`plink --merge-list Small_Bims.txt --out $disease_abbr\_TN_TCGA_All`;

chdir "$RNA_Path";
`mv $RNA_Path/maps/$disease_abbr\_TN_TCGA_All.* $RNA_Path`;

#To save space
`rm -r $RNA_Path/bases`;
`rm -f $RNA_Path/maps/*.bim`;
`rm -f $RNA_Path/maps/*.bed`;
`rm -f $RNA_Path/maps/*.fam`;
`rm -f $RNA_Path/maps/*.log`;

#Split Bed into peds and move them into the dir peds
my @chrs=(1..23);
my @plink_cmds;
while(my $chr=<@chrs>)
{
    push @plink_cmds,"plink --bfile $disease_abbr\_TN_TCGA_All --chr $chr --recode --out $RNA_Path/peds/$chr";
}

mce_map
{
    system("$plink_cmds[$_]");
}0..$#plink_cmds;

`rm -f peds/*.map`;
`rm -f peds/*.log`;

print "All Jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
