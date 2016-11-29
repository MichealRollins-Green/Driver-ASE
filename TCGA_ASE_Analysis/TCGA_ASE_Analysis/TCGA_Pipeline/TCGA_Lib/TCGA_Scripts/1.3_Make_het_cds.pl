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

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $parsing = TCGA_Lib::Parsing_Routines->new;
my $impute_plink = TCGA_Lib::Imputation_Plink->new;
my $Analysispath = realpath("../../Analysis");

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'output|o=s' => \my $imputation,
    'phased_dir|p=s' => \my $phased,
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("1.3");

if($help)
{
    $parsing->usage("1.3");
}

if(!defined $disease_abbr)
{
    print "disease type was not entered!\n";
    $parsing->usage("1.3");
}

my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";

if (!(-d $RNA_Path))
{
    print "$RNA_Path does not exist. Either it was deleted, moved or the scripts required for it have not been run.\n";
    exit;
}

if(!defined $imputation)
{
    $imputation = "$RNA_Path/Impute2Plink";
    print STDERR "using $imputation as default\n"; 
}

unless($phased = "$RNA_Path/phased")
{
    $phased = realpath("../$phased");
}

my $Impute2out="$Analysispath/$disease_abbr/phased_imputed_raw_out";

chdir "$RNA_Path";

`mkdir -p $imputation` unless(-d "$imputation");
`rm -f $RNA_Path/$imputation/*`;

unless("$RNA_Path/Impute2Plink")
{
    $imputation = realpath("../$imputation");
}

my @imputeout = $parsing->get_only_files_in_dir("$Impute2out");
@imputeout = grep {/haps/}@imputeout;
open(I2O,">$RNA_Path/imputed_haps") or die("Can't open file for output:$!");
for(my $i = 0;$i < scalar(@imputeout);$i++)
{
    print I2O "$imputeout[$i]\n";
}

#impute_2_out(path to phased_imputed_raw_out,path to RNA_Seq_Analysis directory,file that contains a list of haps files from directory phased_imputed_raw_out,user defened directory from command line or default directory from script 1.2,user defened directory from command line or default directory)
my @impute_out_cmds = $impute_plink->impute_2_out("$Impute2out","$RNA_Path","$RNA_Path/imputed_haps","$phased","$imputation");

mce_map
{
    system("$impute_out_cmds[$_]");
}0..$#impute_out_cmds;

`rm $RNA_Path/Impute2Plink/*.sample`;
my @merge_list = $parsing->get_only_files_in_dir("$imputation");
@merge_list = grep {/bim/}@merge_list;
open(ML,">$imputation/Merged_list.txt") or die("Can't open file for output:$!");
for(my $i=0;$i<scalar(@merge_list);$i++)
{
 my $bim_file = $merge_list[$i];   
 $bim_file =~ s/.bim//;
 print ML "$bim_file\n";
}
close(ML);

mce_map_f
{
    chomp($_);
    #Assign_ChrPos_MissingId_SNPs(user defened directory from command line or default directory/bim file)
    $impute_plink->Assign_ChrPos_MissingId_SNPs("$imputation/$_.bim");
}"$imputation/Merged_list.txt";

chdir "$imputation";
print "Now doing plink on Merged_list.txt\n";
`plink --threads 3 --merge-list Merged_list.txt --out $disease_abbr\_TN_TCGA_Imputation`;

chdir "$RNA_Path";
`mv $imputation/$disease_abbr\_TN_TCGA_Imputation.* $RNA_Path`;
`rm -r $imputation`;

mkdir "$RNA_Path/cds_plink" unless(-d "$RNA_Path/cds_plink");
`rm -f $RNA_Path/cds_plink/*`;
chdir "$RNA_Path/cds_plink";

print "Now running Extract_Ind_Het_Genos\n";
my @chrs = (1..23);
my @Het_cmds;
while(my $chr=<@chrs>)
{
    `plink --bfile $RNA_Path/$disease_abbr\_TN_TCGA_Imputation --chr $chr --make-bed --out $chr`;
    #Extract_Ind_Het_Genos(chr#,chr file made in the cds_plink directory,21.phased.sample from the user defened directory from command line or default directory from script 1.2,cds_plink directory)
    push @Het_cmds, $impute_plink->Extract_Ind_Het_Genos("chr$chr","$RNA_Path/cds_plink/$chr","$phased/21.phased.sample","$RNA_Path/cds_plink");
};
mce_map
{
    system("$Het_cmds[$_]");
}0..$#Het_cmds;

print "Now running Append_snplist_chrs\n";
mce_map_f
{
    chomp($_);
    if($_=~/TCGA/)
    {
        #Get rid of the first two lines;
        my @as=split("\t| ",$_);
        my $ID=$as[1];
        #Append_snplist_chrs parameters: $ID $snplist_dir and $outdir;
        #Output file name will be $ID.snplist;
        #Append_snplist_chrs(TCGA ID,cds_plink directory)
        $impute_plink->Append_snplist_chrs("$ID","$RNA_Path/cds_plink");
    }
}"$phased/21.phased.sample";

print "Now running merge_tcga_snplist\n";
#merge_tcga_snplist(path to RNA_Seq_Analysis directory)
$impute_plink->merge_tcga_snplist("$RNA_Path");

print "Now running Get_Uniq_SNPS_From_SNPlists\n";
#ALL_SNPLIST.txt is not an actual file instead it is a tag that is used in this subroutine to get files with ALL_SNPLIST.txt in it
#There should be not path included int as it will look for the path in the file name
#Get_Uniq_SNPS_From_SNPlists(cds_plink directory,file tag,output file)
$impute_plink->Get_Uniq_SNPS_From_SNPlists("$RNA_Path/cds_plink","ALL_SNPLIST.txt","all_snplist.txt");

print "Now putting all snp lists into file all_chrs_snplist.txt\n";
$parsing->vlookup("$RNA_Path/cds_plink/all_snplist.txt",1,"$RNA_Path/$disease_abbr\_TN_TCGA_Imputation.bim",2,1,"y","$RNA_Path/all_chrs_snplist.txt");

print "Now running Snplist_Query_Haps_for_Beds\n";
mce_map
{
    my $chr=$_;
    #Snplist_Query_Haps_for_Beds(file that conains all the chr snplists,21.phased.sample in the user defened directory from command line or default directory from script 1.2,path to phased_imputed_raw_out,output bed tag,chr#)
    $impute_plink->Snplist_Query_Haps_for_Beds("$RNA_Path/all_chrs_snplist.txt","$phased/21.phased.sample","$Impute2out","$RNA_Path/cds_plink/All_Hets_Haps","$chr"); 
 }1..23;

print "Now running Split_Bed_Haps_for_individuals\n";
#Split_Bed_Haps_for_individuals(cds_plink directory)
my @split_bed_haps_cmds = $impute_plink->Split_Bed_Haps_for_individuals("$RNA_Path/cds_plink");

mce_map
{
    system("$split_bed_haps_cmds[$_]");
}0..$#split_bed_haps_cmds;

chdir "$RNA_Path/cds_plink";

`rm -f *chr*.snplist *chr*.log *chr*.txt`;
`rm -f *.bim *.ped *.fam`;
`rm -f *chr*.o *chr*.ss *chr*.e`;

#remove non-sorted bed
`ls |grep Het|grep sorted -v|xargs rm`;

chdir $RNA_Path;

mkdir "$RNA_Path/cds_sorted";
`rm -f $RNA_Path/cds_sorted/*`;

print "Now running Split_BedHaps\n";
mce_map
{
    my $chr = $_;
    #Split_BedHaps will only keep Het SNPs;
    #Due to the large number of opened files (sample_size*23=1061*23)
    #Linux system reject the following procedure!
    #Thus, we only open no more than 400 file for writting!
    #Split_BedHaps(All_Hets_Haps.chr#.sorted.bed file in the cds_plink directory,cds_sorted directory,chr#,filehandle#)
    $impute_plink->Split_BedHaps("$RNA_Path/cds_plink/All_Hets_Haps.chr$chr.sorted.bed","$RNA_Path/cds_sorted","$chr","400");
}1..23;

chdir "$RNA_Path/cds_sorted";

print "Now running Merge_All_Chrs_4_Single_TCGA_Bed\n";
#Merge_All_Chrs_4_Single_TCGA_Bed(21.phased.sample file in the the user defened directory from command line or default directory from script 1.2,cds_plink directory,cds_sorted directory)
$impute_plink->Merge_All_Chrs_4_Single_TCGA_Bed("$phased/21.phased.sample","$RNA_Path/cds_plink","$RNA_Path/cds_sorted");

chdir "$RNA_Path";
`rm -f cds_sorted/*chr*`;

`rm -f cds_sorted/*.e cds_sorted/*.o cds_sorted/*.ss`;

mkdir "cds_sorted_tmp";
`rm -f cds_sorted_tmp/*`;

my @cds_sort_files = `ls cds_sorted`;

print "Now running Change_CD_ChrInfo_4_Mpileup\n";
mce_map
{
    chomp(my $TCGA=$_);
    my $rd=rand();
    $rd=substr($rd,2,6);
    #Change_CD_ChrInfo_for_Mpileup(bed file,cds_sorted directory,cds_sorted_tmp/bed file)
    $impute_plink->Change_CD_ChrInfo_for_Mpileup("$TCGA","$RNA_Path/cds_sorted","$RNA_Path/cds_sorted_tmp/$TCGA");
}@cds_sort_files;

#If there are some errors in cds_sorted or cds_sorted_tmp
#Use info in the tar.gz to get right cd files!
chdir "$RNA_Path/cds_plink";
mkdir "$Analysispath/$disease_abbr/$disease_abbr\_finished_analysis_RNA";
print "Now compressing All_Het_Haps and snplists files from cds_plink directory.\n";
`tar -zcvf $disease_abbr\_Imputation_Haps.tar.gz All_Hets_Haps.chr*.sorted.bed TCGA*.snplist`;
`mv $disease_abbr\_Imputation_Haps.tar.gz $Analysispath/$disease_abbr/$disease_abbr\_finished_analysis_RNA`;
chdir $RNA_Path;
`rm -rf cds_plink`;

print "All Jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;