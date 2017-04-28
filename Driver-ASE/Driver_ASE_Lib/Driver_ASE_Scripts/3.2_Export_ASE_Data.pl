#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$Bin/..";
use Gene_Level_ASE_Analysis;
use Parsing_Routines;
use Cwd 'realpath';
use Cwd;
use MCE::Map;
use Getopt::Long;
use strict;

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $ase_analysis = Driver_ASE_Lib::Gene_Level_ASE_Analysis->new;
my $parsing = Driver_ASE_Lib::Parsing_Routines->new;

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage;

if($help)
{
    $parsing->usage;
}

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $tables = "$disease_abbr\_tables";
my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";
my $ase = "$RNA_Path/ase";
my $gene_level = "gene_level";
my ($bad_cnvs,$bad_snps_bed) = ("bad_cnvs","bad_snps_bed");
my $matrix = "matrix";
my ($problem_snps,$problem_cnvs) = ("problem_snps","problem_cnvs");
my $finished_RNA = "$disease_abbr\_finished_analysis_RNA";
my ($temp,$cds_sorted_ase,$logs) = ("temp","cds_sorted_ase","logs");

if(!defined $disease_abbr)
{
    print "disease was not entered!\n";
    $parsing->usage;
}

#Checks if there is no Database directory
if(!(-d "$database_path"))
{
    print STDERR "$database_path does not exist, it was either moved, renamed, deleted or has not been downloaded.\nPlease check the README.md file on the github page to find out where to get the Database directory.\n";
    exit;
}

#Checks if there is no Analysis directory
if (!(-d "$Analysispath"))
{
    print STDERR "$Analysispath does not exist, it was either deleted, moved or renamed.\n";
    print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
    exit;
}
elsif(!(-d "$Analysispath/$disease_abbr"))
{
    print STDERR "$Analysispath/$disease_abbr does not exist, it was either deleted, moved or renamed.\n";
    print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
    exit;
}

if(-d "$Analysispath/$disease_abbr/$tables")
{
   if (!(-e "$Analysispath/$disease_abbr/$tables/final_downloadtable_$disease_abbr\_RNA-Seq.txt"))
    {
        print STDERR "$Analysispath/$disease_abbr/$tables/final_downloadtable_$disease_abbr\_RNA-Seq.txt does not exist. It was either moved, renamed or deleted,\n";
        print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
        exit;
    } 
}
else
{
    print STDERR "$Analysispath/$disease_abbr/$tables does not exist. It was either moved, renamed or deleted,\n";
    print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
    exit;
}

if(!(-d $RNA_Path))
{
    print STDERR "$RNA_Path does not exist. Either it was deleted, moved or renamed.\n";
    print STDERR "Please run script 1.0_Prep_SNPs_for_Imputation_and_Plink.pl.\n";
    exit;
}

if(!(-d "$ase"))
{
    print STDERR "The directory $ase does not exist. It was moved renamed or deleted.\n";
    print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
    exit;
}

chdir $RNA_Path;

if (!(-d "$ase/$gene_level"))
{
    print STDERR "The directory $ase/$gene_level does not exist. It was moved, renamed or deleted.\n";
    print STDERR "Please run script 3.1_ASE_Analysis.pl.\n";
    exit;
}

`ls $ase/$gene_level/ > $ase/gene_level_ases.txt`;
$parsing->vlookup("$ase/gene_level_ases.txt",1,"$ase/RNA_seq_id_lookup_snp6_bed.txt",1,2,"y","$ase/All_ASEs_grep_non_NaN.txt");
`grep NaN -v $ase/All_ASEs_grep_non_NaN.txt > $ase/All_ASEs.txt`;

open(ASE,"$ase/All_ASEs.txt") or die "Can't open $ase/All_ASEs.txt: $!\n";

my @ase = <ASE>;

close(ASE);

open(ASE,">$ase/All_ASEs.txt") or die "Can't open $ase/All_ASEs.txt: $!\n";

foreach my $uuid(@ase)
{
    chomp($uuid);
    my $get_id = [split("\t",$uuid)]->[0];
    $get_id = [split("\\.",$get_id)]->[0];
    
    print ASE $uuid,"\t",$get_id,"\n";
}
close(ASE);

$parsing->vlookup("$ase/All_ASEs.txt",3,"$Analysispath/$disease_abbr/$tables/final_downloadtable_$disease_abbr\_RNA-Seq.txt",1,"6,3","y","$ase/ASEs4lookup_pull_column.txt");

$parsing->pull_column("$ase/ASEs4lookup_pull_column.txt","1,2,4,5","$ase/ASEs4lookup_sort.txt");
`sort -k 2,2 -k 4,4n $ase/ASEs4lookup_sort.txt > $ase/ASEs4lookup_with_bad_CNVs.txt`;

#pull_TN acts as a SAS lag function or retain function
#pull_TN(ase bad cnv file, bad_cnvs directory,bad_snps_bed directory,output file)
$ase_analysis->pull_TN("$ase/ASEs4lookup_with_bad_CNVs.txt","$RNA_Path/$bad_cnvs","$RNA_Path/$bad_snps_bed","$RNA_Path/have_cnv_snps.txt",$disease_abbr);

#functionality of have_cnv_snp is to keep IDs having both bad_cnv and bad_snps
#have_cnv_snp(output file created in pull_TN,output file)
$ase_analysis->have_cnv_snp("$RNA_Path/have_cnv_snps.txt","$RNA_Path/have_cnv_snps_sort.txt");

`sort -k 2,2 $RNA_Path/have_cnv_snps_sort.txt > $RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt`;

#Now link Tumor and Normals gene_level files with TCGA IDs (12 chars)
#Output files will be tumor_ff and normal_ff
#mk_files(the TN snp cnv lookup file created from the sort command, path to the user defined directory from the command line in script 3.0 or default ase)
$ase_analysis->mk_files("$RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt","$ase","$gene_level");

# tumor
#matricize(index file,output file,index column to print,output column,path to RNA_Seq_Analysis directory)
$parsing->matricize("tumor_ff","tumor_ff",4,5,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/tumor.a`;

$parsing->matricize("tumor_ff","tumor_ff",4,6,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/tumor.b`;

$parsing->matricize("tumor_ff","tumor_ff",4,7,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/tumor.p`;

$parsing->matricize("tumor_ff","tumor_ff",4,10,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/tumor.mm`;

$parsing->matricize("tumor_ff","tumor_ff",4,11,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/tumor.agree`;

$parsing->matricize("tumor_ff","tumor_ff",4,12,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/tumor.num_snps`;

`mv $RNA_Path/rowlabels.txt $ase/$matrix`;
`mv $RNA_Path/collabels.txt $ase/$matrix/tumor.collabels`;

# normal
$parsing->matricize("tumor_ff","normal_ff",4,5,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/normal.a`;

$parsing->matricize("tumor_ff","normal_ff",4,6,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/normal.b`;

$parsing->matricize("tumor_ff","normal_ff",4,7,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/normal.p`;

$parsing->matricize("tumor_ff","normal_ff",4,10,"$RNA_Path");
#mm is the number of reads (A,C,G,T)
`mv $RNA_Path/matrix.tab $ase/$matrix/normal.mm`;

$parsing->matricize("tumor_ff","normal_ff",4,11,"$RNA_Path");
#agree: proportion of SNPs with 5+ counts in agreement on direction
#and number of informative snps (5+ counts)
`mv $RNA_Path/matrix.tab $ase/$matrix/normal.agree`;

$parsing->matricize("tumor_ff","normal_ff",4,12,"$RNA_Path");
`mv $RNA_Path/matrix.tab $ase/$matrix/normal.num_snps`;
`mv $RNA_Path/collabels.txt $ase/$matrix/normal.collabels`;
 
 #problem zones -> pull bed from rowlabels and make lists of genes
$parsing->vlookup("$ase/$matrix/rowlabels.txt",1,"$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",4,"1,2,3,4,5,6,7,8,9,10,11,12","n","$ase/$matrix/rowlabels.bed");

#rowlabels.bed will be used by script run_problem_cnvs and run_problem_snps!
#run_problem_cnvs(bad_cnvs directory,the TN snp cnv lookup file created from the sort command,path to RNA_Seq_Analysis directory,matrix directory)
$ase_analysis->run_problem_cnvs("$RNA_Path/$bad_cnvs","$RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt","$ase","$matrix","$problem_snps");
#run_problem_snps(bad_snps_bed directory,the TN snp cnv lookup file created from the sort command,path to RNA_Seq_Analysis directory,matrix directory)
$ase_analysis->run_problem_snps("$RNA_Path/$bad_snps_bed","$RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt","$ase","$matrix","$problem_cnvs");

open(TUM,"tumor_ff") or die "Can't open tumor_ff: $!\n";
open(CNV,">cnv_ff") or die "Can't open cnv_ff:$!\n";

#changes gene_level in path name to problem_cnvs
while (my $r = <TUM>)
{
    chomp($r);
    my @a = split("\t",$r);
    $a[0] =~ s/$gene_level/$problem_cnvs/;
    print CNV join("\t",@a), "\n";
}
close(TUM);
close(CNV);

#matricize_version_two(index_file,print_file,index_column,print_index_column,print_output_column)
$parsing->matricize_version_two("tumor_ff","cnv_ff",4,1,2);

`mv matrix.tab $ase/$matrix/problem.cnv`;

open(TUM,"tumor_ff") or die "Can't open tumor_ff: $!\n";
open(SNP,">snp_ff") or die "Can't open cnv_ff: $!\n";

#changes gene_level in path name to problem_snps
while (my $r = <TUM>)
{
    chomp($r);
    my @a = split("\t",$r); 
    $a[0] =~ s/$gene_level/$problem_snps/;
    print SNP join("\t",@a), "\n";
}
close(SNP);
close(TUM);

$parsing->matricize_version_two("tumor_ff","snp_ff",4,1,2);

`mv matrix.tab $ase/$matrix/problem.snp`;

my @del_files = $parsing->get_only_files_in_dir("$RNA_Path");

@del_files = grep{!/\.bim/ and !/\.fam/ and !/\.log/ and !/\.bed/}@del_files;

for(my $i = 0;$i < scalar(@del_files);$i++)
{
    `rm "$del_files[$i]"`;
}

chdir $ase;

undef @del_files;
@del_files = $parsing->get_only_files_in_dir("$ase");

for(my $i = 0;$i < scalar(@del_files);$i++)
{
    `rm "$del_files[$i]"`;
}

chdir $RNA_Path;

if (-e "$disease_abbr\_ASE_Rst.tar.gz")
{
    print "Now running tar -zcvkf $disease_abbr\_ASE_Rst.tar.gz --exclude \"$ase/$temp\" --exclude \"$ase/$cds_sorted_ase\" --exclude \"$ase/$logs\" $bad_snps_bed $bad_cnvs $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*.\n";
    
    `tar -zcvkf $disease_abbr\_ASE_Rst.tar.gz --exclude "$ase/$temp" --exclude "$ase/$cds_sorted_ase" --exclude "$ase/$logs" $bad_snps_bed $bad_cnvs $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*`;
}
else
{
    print "Now running tar -zcvf $disease_abbr\_ASE_Rst.tar.gz --exclude \"$ase/$temp\" --exclude \"$ase/$cds_sorted_ase\" --exclude \"$ase/$logs\" $bad_snps_bed $bad_cnvs $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*.\n";
    
    `tar -zcvf $disease_abbr\_ASE_Rst.tar.gz --exclude "$ase/$temp" --exclude "$ase/$cds_sorted_ase" --exclude "$ase/$logs" $bad_snps_bed $bad_cnvs $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*`;
}

mkdir "$Analysispath/$disease_abbr/$finished_RNA" unless(-d "$Analysispath/$disease_abbr/$finished_RNA");

`mv $disease_abbr\_ASE_Rst.tar.gz $Analysispath/$disease_abbr/$finished_RNA`;

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
