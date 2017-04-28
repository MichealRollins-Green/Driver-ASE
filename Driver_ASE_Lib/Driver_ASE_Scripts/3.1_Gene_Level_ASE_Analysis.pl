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

if(!defined $disease_abbr)
{
    print "disease type was not entered!\n";
    $parsing->usage;
}

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";
my $ase = "$RNA_Path/ase";
my $mpileups_path = "$ase/rna_mpileups";
my ($cds_sorted_ase,$matrix,$logs,$ase_counts) = ("cds_sorted_ase","matrix","logs","ase_counts");
my $cds_sorted = "cds_sorted";
my $gene_level = "gene_level";

#Checks if there is no Database directory
if(!(-d "$database_path"))
{
    print STDERR "$database_path does not exist, it was either moved, renamed, deleted or has not been downloaded.\nPlease check the README.md file on the github page to find out where to get the Database directory.\n";
    exit;
}

#Checks if there is no Analysis directory
if(!(-d "$Analysispath"))
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

if (!(-d $RNA_Path))
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

chdir "$ase";

mkdir "$cds_sorted_ase" unless(-d "$cds_sorted_ase");
mkdir "$matrix" unless(-d "$matrix");
mkdir "$logs" unless(-d "$logs");
mkdir "$ase_counts" unless(-d "$ase_counts");

`ls $ase/$cds_sorted_ase > $ase/done_sorted_ase.txt`;

if (!(-d "$RNA_Path/$cds_sorted"))
{
    print STDERR "The directory $RNA_Path/$cds_sorted does not exist. It was moved, renamed or deleted\n";
    print STDERR "Please run script 1.3_Make_het_cds.pl.\n";
    exit;
}

`ls $RNA_Path/$cds_sorted > $ase/cds_for_change_beds.txt`;

#checks for bed files that have not yet been done and processes them for cds_sorted_ase. This only accounts for beds that are not in the directory and any incomplete beds will not be processed so if there are any beds suspected to not be complete, delete them before running this script.
$parsing->vlookup("$ase/cds_for_change_beds.txt",1,"$ase/done_sorted_ase.txt",1,1,"y","cds_for_beds_gl_grep.txt");
`grep NaN cds_for_beds_gl_grep.txt > cds_for_beds_gl.txt`;

if(-s "cds_for_beds_gl.txt" > 0)
{
    #removes rsIDs in cds_beds for later analysis of gene level ASE!
    #Change all beds in cds_sorted for later calling of gene level ase!
    #change_beds_gl_ase(cds file with list of beds from the cds_sorted directory,path to the RNA_Seq_Analysis directory,user defined directory from the command line in script 3.0 or default ase)
    $ase_analysis->change_beds_gl_ase("$ase/cds_for_beds_gl.txt","$RNA_Path/$cds_sorted","$ase/$cds_sorted_ase");
}

my @mpileups = `ls $mpileups_path|grep '.'`;
@mpileups = grep{chomp;-s "$mpileups_path/$_";}@mpileups;

my @cds = `ls $RNA_Path/$cds_sorted`;

open(MP,">RNA_seq_id_lookup_snp6_bed.txt") or die "cant open RNA_seq_id_lookup_snp6_bed.txt: $!\n";

foreach my $mp(@mpileups)
{
    my $tcga = [split("\\.",$mp)]->[-1];
    my @bed = grep{/$tcga/}@cds;
    print MP $mp, "\t", $bed[0];
}
close(MP);

#Get already done files in ase_counts and lookup it with RNA_seq4ASEPipeup.txt
#Only run these files absent in the ase_counts dir
`ls $ase_counts > already_done.txt`;
`ls $mpileups_path/ > ase_wish_list.txt`;

$parsing->vlookup("$ase/ase_wish_list.txt",1,"$ase/already_done.txt",1,1,"y","left_ase_grep_NaN.txt");
`grep NaN left_ase_grep_NaN.txt > left_ase_pull_column.txt`;
$parsing->pull_column("left_ase_pull_column.txt",1,"left_ase_counts.txt");

$parsing->vlookup("$ase/RNA_seq_id_lookup_snp6_bed.txt",1,"left_ase_counts.txt",1,1,"y","RNA_seq_id_grep_non_NaN");
`grep NaN -v RNA_seq_id_grep_non_NaN > RNA_seq_id_pull_column.txt`;
$parsing->pull_column("RNA_seq_id_pull_column.txt","1,2","RNA_seq_id_ase_no_tum_norm.txt");

# ase @ snps - also make lookup (tumor_norm_look) for matricization
# compile_ase_no_tum_norm will use
# subroutne pileup_at_cd
# It will import contents of cds_bed files first
#    0: chr1      
#    1: 245042304 
#    2: 245042305 
#    3: T|C       
# Then it will count the two alleles in cds_bed in the sequence line in mpileup output
# The sequence line is in the 6 columns
#    0: chr1                         
#    1: 245344663                    
#    2: N                            
#    3: 29                           
#    4: >>>>>>>>>>>>>>><<>>>><>>>>><<
#    5: FFFJJGJIIJJIIIIFIJBJIIHHHHHHH

# The final results containing the number of ase_counts will be saved in dir ase_counts
# The ase_counts is CNV level count
#compile_ase_no_tum_norm(file that conatins list of mpileups and associated bed files, cds_sorted directory)
$ase_analysis->compile_ase_no_tum_norm("RNA_seq_id_ase_no_tum_norm.txt","$RNA_Path/$cds_sorted","$mpileups_path",$ase_counts);

#get gene level ase_counts

mkdir "$ase/$gene_level" unless(-d "$ase/$gene_level");
`ls $ase/$gene_level > already_done_genes.txt`;
`ls $mpileups_path/ > ases.txt`;

$parsing->vlookup("$ase/ases.txt",1,"$ase/already_done_genes.txt",1,1,"y","$ase/left_gene_grep_NaN.txt");
`grep NaN $ase/left_gene_grep_NaN.txt > $ase/left_gene_pull_collumn.txt`;
$parsing->pull_column("$ase/left_gene_pull_collumn.txt",1,"$ase/left_gene_level_ase.txt");

$parsing->vlookup("$ase/RNA_seq_id_lookup_snp6_bed.txt",1,"$ase/left_gene_level_ase.txt",1,1,"y","$ase/RNA_seq_id_lookup_grep_non_NaN.txt");
`grep NaN -v $ase/RNA_seq_id_lookup_grep_non_NaN.txt > $ase/RNA_seq_id_lookup_pull_column.txt`;
$parsing->pull_column("$ase/RNA_seq_id_lookup_pull_column.txt","1,2","$ase/RNA_seq_id_lookup_comp_gene_faster.txt");
#compile_gene_ase_faster will use bed file saved in the directory cds_ase_sorted.
#compile_gene_ase_faster(file that conatins list of ase_counts and associated bed files,cds_sorted_ase directory,path to the refseq.ucsc.ensembl.mrna.hg9.nr.bed file in the Database directory,user defined directory from the command line in script 3.0 or default ase)
$ase_analysis->compile_gene_ase_faster("$ase/RNA_seq_id_lookup_comp_gene_faster.txt","$ase/$cds_sorted_ase","$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed","$ase","$ase_counts","$gene_level");

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
