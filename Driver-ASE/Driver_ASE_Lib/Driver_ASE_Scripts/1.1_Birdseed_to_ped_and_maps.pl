#!/usr/bin/perl -w

use MCE::Map;
use FindBin qw($Bin);
use lib "$Bin/..";
use Parsing_Routines;
use Imputation_Plink;
use Cwd 'realpath';
use Getopt::Long;
use strict;
use autodie;
no warnings 'once';

my $time = localtime;
print "Script started: $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $parsing = Driver_ASE_Lib::Parsing_Routines->new;
my $impute_plink = Driver_ASE_Lib::Imputation_Plink->new;

GetOptions(
    'cancer|c=s' => \my $cancer_type, #e.g. OV
    'plink|p=s' => \my $plink, #path to plink
    'sampletype|s=i' => \my $sampletype, #0(normal) 1(tumor) 2(normal/tumor)
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("1.1");

if ($help)
{
    $parsing->usage("1.1");
}

if (!defined $cancer_type)
{
    print STDERR "Cancer type was not entered!\n";
    $parsing->usage("1.1");
}

if (defined $sampletype)
{
    if ($sampletype != 0 and $sampletype != 1 and $sampletype != 2)
    {
	print STDERR "Invalid input type! The input type must be 0|1|2.";
	$parsing->usage("1.1");
    }
}

if (!defined $plink)
{
    $plink = "plink";
}

$parsing->CheckSoftware("$plink");

my $Driver_ASE_Dir = realpath("../../");
my $Analysispath = realpath("../../Analysis");
my $RNA_Path = "$Analysispath/$cancer_type/RNA_Seq_Analysis";
my $affy_dir = "affy6";
my $bases_dir = "$RNA_Path/bases";
my $SNP6_Raw_Files_Dir = "$Analysispath/$cancer_type/SNP6/Genotypes";
my $map_dir = "maps";
my $ped_dir = "peds";
my $database_path = "$Driver_ASE_Dir/Database";
my $snp6_cd_file = "snp6.cd.txt";

$parsing->check_directory_existence("$database_path","$Analysispath","$Analysispath/$cancer_type","$SNP6_Raw_Files_Dir","$RNA_Path","$RNA_Path/$affy_dir","$RNA_Path/$affy_dir/$snp6_cd_file"); #check if directories or files exist
$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

chdir "$RNA_Path";

`mkdir -p "$bases_dir"` unless(-d "$bases_dir");
`rm -rf $bases_dir/*`;

#gets all of the files in the the directory that was specified to download the SNP-Array data
my @birdseed = $parsing->get_only_files_in_dir("$SNP6_Raw_Files_Dir");

mce_map
{
    chomp ($_);
   #Prepare_Genos(Genotype files,$SNP6_Raw_Files_Dir (directory where Genotype files were downloaded),$RNA_Path/$affy_dir/$snp6_cd_file (snp6_cd.txt file created from 1.0),$bases_dir (path to directory with default name bases),$sampletype(e.g. 0(normal), 1(tumor) or 2(normal/tumor)))
    $impute_plink->Prepare_Genos("$_","$SNP6_Raw_Files_Dir","$RNA_Path/$affy_dir/$snp6_cd_file","$bases_dir",$sampletype);
}@birdseed;

#remove the .t files in the $bases_dir
my @t_files = `ls "$bases_dir"`;
chdir "$bases_dir";
for (my $i = 0;$i < scalar(@t_files);$i++)
{
    if ($t_files[$i] =~ /\.t/)
    {
        chomp($t_files[$i]);
        `rm '$bases_dir'/$t_files[$i]`;  
    }
}

chdir "$RNA_Path";

 my @files = $parsing->get_only_files_in_dir("$bases_dir");
@files = grep{!/^\./ && -f "$bases_dir/$_"}@files;
 my $selected_file = $files[-1];
 #pull_column(file to pull column from,column(s) to pull, output file)
 $parsing->pull_column("$bases_dir/$selected_file","1","$cancer_type\_t");

 mkdir "$map_dir" unless(-d "$RNA_Path/$map_dir");
 `rm -rf $RNA_Path/$map_dir/*`;

 $parsing->vlookup("$cancer_type\_t",1,"$RNA_Path/$affy_dir/$snp6_cd_file",4,"1,3","y","map_for_sort");

 `sort -k 2,2 -k 3,3n map_for_sort > map_for_pull_column`;
 $parsing->pull_column("map_for_pull_column","2,3,1","mk_map_file");
 #mk_map(file used to make maps,$RNA_Path (path to directory where RNA-Seq analysis data is stored))
 $impute_plink->mk_map("mk_map_file","$RNA_Path/$map_dir");

 my @maps = $parsing->get_only_files_in_dir("$RNA_Path/$map_dir");
 @maps = grep{!/^\./ and /^[0-9x]+.map$/i}@maps;
 open (CHR,">chrs");
 for (my $i=0;$i<scalar(@maps);$i++)
{
    print CHR $maps[$i], "\n";
}
 close (CHR);

mkdir "$RNA_Path/$ped_dir" unless(-d "$RNA_Path/$ped_dir");
`rm -rf $RNA_Path/$ped_dir/*`;

#mk_ped_list_4_plink(file with chr.map,$bases_dir (directory with default name of bases),outfile for peds,$RNA_Path (path to directory where RNA-Seq analysis data is stored))
 $impute_plink->mk_ped_list_4_plink("chrs","$bases_dir","ped_list","$RNA_Path");

  mce_map_f
  {
    chomp($_);
    my @lines = split("\t",$_);
    #only keep normal samples for plink
    my $tn = [split("-",$lines[1])]->[-1];
    $tn =~ s/[a-zA-Z]$//;
    #Make_Plink_Bed(chr#,sample(e.g. 10),sample id,path to directory where RNA-Seq analysis data is stored,$bases_dir (directory with default name of bases),$map_dir (path to directory where map files are stored),path/command for plink)
    $impute_plink->Make_Plink_Bed("$lines[0]","$lines[1]","$lines[2]","$lines[3]","$bases_dir","$map_dir","$plink") if $tn > 9;
  }"$RNA_Path/ped_list";

  my @maps_bim = $parsing->get_only_files_in_dir("$RNA_Path/$map_dir");
  @maps_bim = grep {/bim/}@maps_bim;
  open (BIMS,">$RNA_Path/$map_dir/Small_Bims.txt");
  for (my $i=0;$i<scalar(@maps_bim);$i++)
  {
    $maps_bim[$i] =~ s/\.bim//;
    print BIMS "$maps_bim[$i]\n";
  }
  close (BIMS);

 chdir "$RNA_Path/$map_dir";
 #Merge all bims in dir maps with plink
 `$plink --merge-list Small_Bims.txt --out $cancer_type\_TN_TCGA_All`;

chdir "$RNA_Path";
`mv $RNA_Path/$map_dir/$cancer_type\_TN_TCGA_All.* $RNA_Path`;

#To save space
`rm -rf $bases_dir` unless(!(-d "$bases_dir"));
my @del_files = `ls $RNA_Path/$map_dir`;

 for (my $i = 0;$i < scalar(@del_files);$i++)
 {
    if ($del_files[$i] =~ /.bim$/ || $del_files[$i] =~ /.bed$/ || $del_files[$i] =~ /.fam$/ || $del_files[$i] =~ /.log$/ || $del_files[$i] =~ /.ped$/ || $del_files[$i] =~ /.pro$/)
    {
        `rm -f $RNA_Path/$map_dir/$del_files[$i]`;
    }
 }

#Split Bed into peds and move them into the dir peds
my @chrs=(1..23);
my @plink_cmds;
while (my $chr=<@chrs>)
{
    push @plink_cmds,"$plink --bfile $cancer_type\_TN_TCGA_All --chr $chr --recode --out $RNA_Path/$ped_dir/$chr";
}

mce_map
{
    system("$plink_cmds[$_]");
}0..$#plink_cmds;

`rm -f $ped_dir/*.map`;
`rm -f $ped_dir/*.log`;

print "All Jobs have finished for $cancer_type.\n";

$time = localtime;
print "Script finished: $time.\n";

exit;
