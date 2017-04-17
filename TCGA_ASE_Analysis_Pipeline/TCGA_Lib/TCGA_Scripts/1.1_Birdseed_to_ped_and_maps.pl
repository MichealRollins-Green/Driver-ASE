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
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage;

if($help)
{
    $parsing->usage;
}

my $Analysispath = realpath("../../Analysis");
my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";
my $affy_dir = "affy6";
my $bases_dir = "$RNA_Path/bases";
my $SNP6_Raw_Files_Dir = "$Analysispath/$disease_abbr/SNP6/Genotypes";
my $map_dir = "maps";
my $ped_dir = "peds";

if(!defined $disease_abbr)
{
    print "disease type was not entered!\n";
    $parsing->usage;
}

#Checks if there is no Analysis directory
if (!(-d "$Analysispath"))
{
    print STDERR "$Analysispath does not exist. It was either deleted, moved or renamed.\n";
    print STDERR "Please run script 0_Download_SNPArray_From_GDC.pl.\n";
    exit;
}
elsif(!(-d "$Analysispath/$disease_abbr"))
{
    print STDERR "$Analysispath/$disease_abbr does not exist. It was either deleted, moved or renamed.\n";
    print STDERR "Please run script 0_Download_SNPArray_From_GDC.pl.\n";
    exit;
}

if (!(-d $RNA_Path))
{
    print STDERR "$RNA_Path does not exist. Either it was deleted, moved or renamed.\n";
    print STDERR "Please run script 1.0_Prep_SNPs_for_Imputation_and_Plink.pl.\n";
    exit;
}

if(!(-d "$RNA_Path/$affy_dir"))
{
    print STDERR "The directory $RNA_Path/$affy_dir does not exist. It was moved, renamed or deleted.\n";
    print STDERR "Please run script 1.0_Prep_SNPs_for_Imputation_and_Plink.pl.\n";
    exit;
}

chdir "$RNA_Path";

`mkdir -p "$bases_dir"` unless(-d "$bases_dir");
`rm -rf $bases_dir/*`;

if (!(-d "$SNP6_Raw_Files_Dir"))
{
    print STDERR "The directory $SNP6_Raw_Files_Dir does not exist. It was moved, renamed or deleted";
    print STDERR "Please run script 0_Download_SNPArray_From_GDC.pl.\n";
    exit;
}

#gets all of the files in the the directory that was specified to download the SNP-Array data
my @birdseed = $parsing->get_only_files_in_dir("$SNP6_Raw_Files_Dir");

mce_map
{
    chomp ($_);
   #Prepare_Genos(Genotype files,directory where Genotype files were downloaded,snp6.cd.txt file created from 1.0,user or default output directory from command line,sample type(e.g. 0(normal), 1(tumor) or 2(normal/tumor))
    $impute_plink->Prepare_Genos("$_","$SNP6_Raw_Files_Dir","$RNA_Path/$affy_dir/snp6.cd.txt","$bases_dir",2);
}@birdseed;

#remove the .t files in the $bases_dir
my @t_files = `ls "$bases_dir"`;
chdir "$bases_dir";
for(my $i = 0;$i < scalar(@t_files);$i++)
{
    if($t_files[$i] =~ /\.t/)
    {
        chomp($t_files[$i]);
        `rm '$bases_dir'/$t_files[$i]`;  
    }
}

chdir "$RNA_Path";

 my @files = $parsing->get_only_files_in_dir("$bases_dir");
    @files=grep{!/^\./ && -f "$bases_dir/$_"}@files;
 my $selected_file = $files[-1];
 #pull_column(file to pull column from,column(s) to pull, output file)
 $parsing->pull_column("$bases_dir/$selected_file","1","$disease_abbr\_t");

 mkdir "$map_dir" unless(-d "$RNA_Path/$map_dir");
 `rm -rf $RNA_Path/$map_dir/*`;

 $parsing->vlookup("$disease_abbr\_t",1,"$RNA_Path/$affy_dir/snp6.cd.txt",4,"1,3","y","map_for_sort");

 `sort -k 2,2 -k 3,3n map_for_sort > map_for_pull_column`;
 $parsing->pull_column("map_for_pull_column","2,3,1","mk_map_file");
 #mk_map(file used to make maps, path to RNA_Seq_Analysis directory)
 $impute_plink->mk_map("mk_map_file","$RNA_Path/$map_dir");

 my @maps = $parsing->get_only_files_in_dir("$RNA_Path/$map_dir");
 @maps = grep{!/^\./ and /^[0-9x]+.map$/i}@maps;
 open(CHR,">chrs") or die("Can't open filehandle: $!");
 for(my $i=0;$i<scalar(@maps);$i++)
{
     print CHR $maps[$i], "\n";
 }
 close(CHR);

mkdir "$RNA_Path/$ped_dir" unless(-d "$RNA_Path/$ped_dir");
`rm -rf $RNA_Path/$ped_dir/*`;

#mk_ped_list_4_plink(file with chr.map,user or default output directory from command line,outfile for peds,path to RNA_Seq_Analysis directory)
 $impute_plink->mk_ped_list_4_plink("chrs","$bases_dir","ped_list","$RNA_Path");

  mce_map_f
  {
      chomp($_);
      my @lines = split("\t",$_);
      #only keep normal samples for plink
      my $tn = [split("-",$lines[1])]->[-1];
      $tn =~ s/[a-zA-Z]$//;
      #Make_Plink_Bed(chr#,sample(e.g. 10),sample id,path to RNA_Seq_Analysis directory,user or default output directory from command line)
      $impute_plink->Make_Plink_Bed("$lines[0]","$lines[1]","$lines[2]","$lines[3]","$bases_dir","$map_dir") if $tn > 9;
  }"$RNA_Path/ped_list";

  my @maps_bim = $parsing->get_only_files_in_dir("$RNA_Path/$map_dir");
  @maps_bim = grep {/bim/}@maps_bim;
  open(BIMS,">$RNA_Path/$map_dir/Small_Bims.txt") or die("Can't open file for output: $!");
  for(my $i=0;$i<scalar(@maps_bim);$i++)
  {
      $maps_bim[$i] =~ s/\.bim//;
      print BIMS "$maps_bim[$i]\n";
  }
  close(BIMS);

 chdir "$RNA_Path/$map_dir";
 #Merge all bims in dir maps with plink
 `plink --merge-list Small_Bims.txt --out $disease_abbr\_TN_TCGA_All`;

chdir "$RNA_Path";
`mv $RNA_Path/$map_dir/$disease_abbr\_TN_TCGA_All.* $RNA_Path`;

#To save space
`rm -rf $bases_dir` unless(!(-d "$bases_dir"));
my @del_files = `ls $RNA_Path/$map_dir`;

 for(my $i = 0;$i < scalar(@del_files);$i++)
 {
     if ($del_files[$i] =~ /.bim$/)
     {
         `rm -f $RNA_Path/$map_dir/$del_files[$i]`;
     }
     elsif($del_files[$i] =~ /.bed$/)
     {
         `rm -f $RNA_Path/$map_dir/$del_files[$i]`;
     }
     elsif($del_files[$i] =~ /.fam$/)
     {
       `rm -f $RNA_Path/$map_dir/$del_files[$i]`; 
     }
     elsif($del_files[$i] =~ /.log$/)
     {
         `rm -f $RNA_Path/$map_dir/$del_files[$i]`;
     }
     elsif($del_files[$i] =~ /.ped$/)
     {
         `rm -f $RNA_Path/$map_dir/$del_files[$i]`;  
     }
     elsif($del_files[$i] =~ /.pro$/)
     {
         `rm -f $RNA_Path/$map_dir/$del_files[$i]`;  
     }
 }

#Split Bed into peds and move them into the dir peds
my @chrs=(1..23);
my @plink_cmds;
while(my $chr=<@chrs>)
{
    push @plink_cmds,"plink --bfile $disease_abbr\_TN_TCGA_All --chr $chr --recode --out $RNA_Path/$ped_dir/$chr";
}

mce_map
{
    system("$plink_cmds[$_]");
}0..$#plink_cmds;

`rm -f $ped_dir/*.map`;
`rm -f $ped_dir/*.log`;

print "All Jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
