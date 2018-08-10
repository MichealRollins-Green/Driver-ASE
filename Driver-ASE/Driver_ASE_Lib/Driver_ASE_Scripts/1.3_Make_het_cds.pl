#!/usr/bin/perl -w

use MCE::Map;
use FindBin qw($Bin);
use lib "$Bin/..";
use Parsing_Routines;
use Imputation_Plink;
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
my $impute_plink = Driver_ASE_Lib::Imputation_Plink->new;

GetOptions(
    'cancer|c=s' => \my $cancer_type, #e.g. OV
    'plink|p=s' => \my $plink, #path to plink
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("1.3");

if ($help)
{
    $parsing->usage("1.3");
}

if (!defined $cancer_type)
{
    print STDERR "Cancer type was not entered!\n";
    $parsing->usage("1.3");
}

if (!defined $plink)
{
    $plink = "plink";
}

$parsing->CheckSoftware("$plink");

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $RNA_Path = "$Analysispath/$cancer_type/RNA_Seq_Analysis";
my $imputation = "$RNA_Path/Impute2Plink";
my $phased = "$RNA_Path/phased";
my $Impute2out = "$Analysispath/$cancer_type/phased_imputed_raw_out";
my $cds_plink = "cds_plink";
my $cds_sorted = "cds_sorted";
my $cds_tmp = "cds_sorted_tmp";
my $finished_RNA = "$cancer_type\_finished_analysis_RNA";
my $phased_sample = "21.phased.sample";

$parsing->check_directory_existence("$database_path","$Analysispath","$Analysispath/$cancer_type","$Impute2out","$RNA_Path","$phased"); #check if directories or files exist
$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

chdir "$RNA_Path";

`mkdir -p $imputation` unless(-d "$imputation");
`find $RNA_Path/$imputation/* 2>/dev/null |xargs rm -rf`;
#`rm -f $RNA_Path/$imputation/*`;

my @imputeout = $parsing->get_only_files_in_dir("$Impute2out");
@imputeout = grep {/haps/}@imputeout;
open (I2O,">$RNA_Path/imputed_haps");
for (my $i = 0;$i < scalar(@imputeout);$i++)
{
    print I2O "$imputeout[$i]\n";
}
close (I2O);
#impute_2_out($Impute2out (path to directory where raw imputed data was stored in script 1.2),$RNA_Path (path to directory where RNA-Seq Analysis is stored),$RNA_Path/imputed_haps (file that contains a list of haps files from directory where raw imputed is stored),$phased (path to directory with default name phased),$imputation (path to directory where impluted plink data will be stored),$plink (path/command for plink),$phased_sample (phased sample file))
my @impute_out_cmds = $impute_plink->impute_2_out("$Impute2out","$RNA_Path","$RNA_Path/imputed_haps","$phased","$imputation","$plink","$phased_sample");

mce_map
{
    system("$impute_out_cmds[$_]");
}0..$#impute_out_cmds;

`rm $imputation/*.sample`;
my @merge_list = $parsing->get_only_files_in_dir("$imputation");
@merge_list = grep {/bim/}@merge_list;
open (ML,">$imputation/Merged_list.txt");
for (my $i=0;$i<scalar(@merge_list);$i++)
{
 my $bim_file = $merge_list[$i];   
 $bim_file =~ s/.bim//;
 print ML "$bim_file\n";
}
close (ML);

mce_map_f
{
    chomp($_);
    #Assign_ChrPos_MissingId_SNPs($imputation/$_.bim (path to bim files created in impute_2_out))
    $impute_plink->Assign_ChrPos_MissingId_SNPs("$imputation/$_.bim");
}"$imputation/Merged_list.txt";

chdir "$imputation";
print "Now doing plink on Merged_list.txt.\n";
`$plink --threads 3 --merge-list Merged_list.txt --out $cancer_type\_TN_TCGA_Imputation`;

chdir "$RNA_Path";
`mv $imputation/$cancer_type\_TN_TCGA_Imputation.* $RNA_Path`;

mkdir "$RNA_Path/$cds_plink" unless(-d "$RNA_Path/$cds_plink");
`find $RNA_Path/$cds_plink/* 2>/dev/null |xargs rm -rf`;
#`rm -f $RNA_Path/$cds_plink/*`;
chdir "$RNA_Path/$cds_plink";

print "Now running Extract_Ind_Het_Genos.\n";
my @chrs = (1..23);
my @Het_cmds;
while (my $chr=<@chrs>)
{
    `$plink --bfile $RNA_Path/$cancer_type\_TN_TCGA_Imputation --chr $chr --make-bed --out $chr`;
    #Extract_Ind_Het_Genos(chr$chr (chr#),$RNA_Path/$cds_plink/$chr (chr file made in the directory with default name cds_plink),$phased/$phased_sample (path to phased sample file),$RNA_Path/$cds_plink (path to directory with default name cds_plink),$plink (path/command for plink))
    push @Het_cmds, $impute_plink->Extract_Ind_Het_Genos("chr$chr","$RNA_Path/$cds_plink/$chr","$phased/$phased_sample","$RNA_Path/$cds_plink","$plink");
};
mce_map
{
    system("$Het_cmds[$_]");
}0..$#Het_cmds;

print "Now running Append_snplist_chrs.\n";
mce_map_f
{
    chomp($_);
    if ($_=~/TCGA/)
    {
        #Get rid of the first two lines;
        my @as=split("\t| ",$_);
        my $ID=$as[1];
        #Append_snplist_chrs parameters: $ID $snplist_dir and $outdir;
        #Output file name will be $ID.snplist;
        #Append_snplist_chrs($ID (TCGA ID),$RNA_Path/$cds_plink (path to directory with default name cds_plink))
        $impute_plink->Append_snplist_chrs("$ID","$RNA_Path/$cds_plink");
    }
}"$phased/$phased_sample";

print "Now running merge_tcga_snplist.\n";
#merge_tcga_snplist($RNA_Path/$cds_plink (path to directory with default name cds_plink))
$impute_plink->merge_tcga_snplist("$RNA_Path/$cds_plink");

print "Now running Get_Uniq_SNPS_From_SNPlists.\n";
#ALL_SNPLIST.txt is not an actual file instead it is a tag that is used in this subroutine to get files with ALL_SNPLIST.txt in it
#There should be not path included int as it will look for the path in the file name
#Get_Uniq_SNPS_From_SNPlists($RNA_Path/$cds_plink (path to directory with default name cds_plink),file tag,output file)
$impute_plink->Get_Uniq_SNPS_From_SNPlists("$RNA_Path/$cds_plink","ALL_SNPLIST.txt","all_snplist.txt");

print "Now putting all snp lists into file all_chrs_snplist.txt.\n";
$parsing->vlookup("$RNA_Path/$cds_plink/all_snplist.txt",1,"$RNA_Path/$cancer_type\_TN_TCGA_Imputation.bim",2,1,"y","$RNA_Path/all_chrs_snplist.txt");

print "Now running Snplist_Query_Haps_for_Beds.\n";
mce_map
{
    my $chr=$_;
    #Snplist_Query_Haps_for_Beds(file that conains all the chr snplists,$phased/$phased_sample (path to phased sample file),$Impute2out (path to directory where raw imputed data was stored in script 1.2),output bed tag,chr#)
    $impute_plink->Snplist_Query_Haps_for_Beds("$RNA_Path/all_chrs_snplist.txt","$phased/$phased_sample","$Impute2out","$RNA_Path/$cds_plink/All_Hets_Haps","$chr"); 
 }1..23;

print "Now running Split_Bed_Haps_for_individuals.\n";
#Split_Bed_Haps_for_individuals($RNA_Path/$cds_plink (path to directory with default name cds_plink))
my @split_bed_haps_cmds = $impute_plink->Split_Bed_Haps_for_individuals("$RNA_Path/$cds_plink");

mce_map
{
    `$split_bed_haps_cmds[$_]`;
}0..$#split_bed_haps_cmds;

chdir "$RNA_Path/$cds_plink";

system("find $RNA_Path/$cds_plink/ -type f -name \"*chr*.snplist\"-o -name \"*chr*.log\" -o -name \"*chr*.txt\" -o -name \"*.bim\" -o -name \"*.ped\" -o -name \"*.fam\" -o -name \"*chr*.o\" -o -name \"*chr*.ss\" -o -name \"*chr*.e\" 2>/dev/null|xargs rm -rf");

#`rm -f *chr*.snplist *chr*.log *chr*.txt`;
#`rm -f *.bim *.ped *.fam`;
#`rm -f *chr*.o *chr*.ss *chr*.e`;

#remove non-sorted bed
`ls |grep Het|grep sorted -v|xargs rm`;

chdir $RNA_Path;

mkdir "$RNA_Path/$cds_sorted" unless(-d "$RNA_Path/$cds_sorted");
`find $RNA_Path/$cds_sorted/* 2>/dev/null |xargs rm -rf`;
#`rm -f $RNA_Path/$cds_sorted/*`;

print "Now running Split_BedHaps.\n";
mce_map
{
    my $chr = $_;
    #Split_BedHaps will only keep Het SNPs;
    #Due to the large number of opened files (sample_size*23=1061*23)
    #Linux system reject the following procedure!
    #Thus, we only open no more than 400 file for writting!
    #Split_BedHaps($RNA_Path/$cds_plink/All_Hets_Haps.chr$chr.sorted.bed (All_Hets_Haps.chr#.sorted.bed file in the directory with default name cds_plink),$RNA_Path/$cds_sorted (path to directory with default name cds_sorted),chr#,filehandle#)
    $impute_plink->Split_BedHaps("$RNA_Path/$cds_plink/All_Hets_Haps.chr$chr.sorted.bed","$RNA_Path/$cds_sorted","$chr","400");
}1..23;

chdir "$RNA_Path/$cds_sorted";

print "Now running Merge_All_Chrs_4_Single_TCGA_Bed.\n";
#Merge_All_Chrs_4_Single_TCGA_Bed($phased/$phased_sample (path to phased sample file created in script 1.2),$RNA_Path/$cds_plink (path to directory with default name cds_plink),$RNA_Path/$cds_sorted (path to directory with default name cds_sorted))
$impute_plink->Merge_All_Chrs_4_Single_TCGA_Bed("$phased/$phased_sample","$RNA_Path/$cds_plink","$RNA_Path/$cds_sorted");

system("find $RNA_Path/$cds_sorted/ -type f -name \"*chr*\" -o -name \"*e\" -o -name \"*o\" -o -name \"*ss\" 2>/dev/null|xargs rm -rf");

chdir "$RNA_Path";
#`rm -f $cds_sorted/*chr*`;

#`rm -f $cds_sorted/*.e $cds_sorted/*.o $cds_sorted/*.ss`;

mkdir "$cds_tmp" unless(-d "$RNA_Path/$cds_tmp");
`find $RNA_Path/$cds_tmp/* 2>/dev/null |xargs rm -rf`;
#`rm -f $cds_tmp/*`;

my @cds_sort_files = `ls $cds_sorted`;

print "Now running Change_CD_ChrInfo_4_Mpileup.\n";
mce_map
{
    chomp(my $TCGA=$_);
    my $rd=rand();
    $rd=substr($rd,2,6);
    #Change_CD_ChrInfo_for_Mpileup($TCGA (bed file),$RNA_Path/$cds_sorted (path to directory with default name cds_sorted),$RNA_Path/$cds_tmp/$TCGA ( path to bed file created in directory with default name cds_sorted_tmp))
    $impute_plink->Change_CD_ChrInfo_for_Mpileup("$TCGA","$RNA_Path/$cds_sorted","$RNA_Path/$cds_tmp/$TCGA");
}@cds_sort_files;

print "All Jobs have finished for $cancer_type.\n";

$time = localtime;
print "Script finished: $time.\n";

exit;
