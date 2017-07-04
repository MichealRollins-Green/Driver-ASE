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
my ($matrix,$logs,$ase_counts) = ("matrix","logs","ase_counts");
my $cds_sorted = "cds_sorted";
my $gene_level = "gene_level";

#Checks if there is no Database directory
if(!(-d "$database_path"))
{
    print STDERR "$database_path does not exist, it was either moved, renamed, deleted or has not been downloaded.\nPlease check the README.md file on the github page to find out where to get the Database directory.\n";
    exit;
}

#Check if the cancer type entered exists with in the file.
open(my $can,"$database_path/Cancer_Types.txt") or die "Can't open Cancer_Types.txt for input: $!\n";
my @can_types = <$can>;
my $line_num = 1;
my $num_of_ctypes = scalar(@can_types);
my $no_count = 0;

print "Locating $disease_abbr...\n";

foreach my $line (@can_types)
{
    chomp($line);
    if ($disease_abbr eq $line)
    {
	print "Found $disease_abbr on line $line_num.\n\nContinuing program.\n\n";
	last;
    }
    else
    {
	print "No $disease_abbr on line $line_num.\n";
	print "Line $line_num was $line.\n\n";
	$no_count += 1;
    }
    $line_num += 1;
}
close ($can);

if ($no_count == $num_of_ctypes)
{
    print "$disease_abbr is not in the Cancer_Types.txt file. Maybe it was misspelled or it does not exits within the file.\n";
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

if (!(-d "$RNA_Path/$cds_sorted"))
{
    print STDERR "The directory $RNA_Path/$cds_sorted does not exist. It was moved, renamed or deleted\n";
    print STDERR "Please run script 1.3_Make_het_cds.pl.\n";
    exit;
}

mkdir "$matrix" unless(-d "$matrix");
mkdir "$logs" unless(-d "$logs");
mkdir "$ase_counts" unless(-d "$ase_counts");

#Get mpileup files;
opendir (MP, $mpileups_path) or die "can not open the directory $mpileups_path: $!\n";
my @mfs = grep{!/^\./ && -f "$mpileups_path/$_"}readdir(MP);
closedir(MP);

#Get cds files;
opendir (CD, "$RNA_Path/$cds_sorted") or die "can not open the directory $RNA_Path/$cds_sorted: $!\n";
my @cds= grep{!/^\./ && -f "$RNA_Path/$cds_sorted/$_"} readdir(CD);
closedir(CD);
@cds = map{
            s/\.bed//;$_;
            }@cds;

#match cds with its corresponding mpileup file;
open(my $T, ">RNA_CDS_mpileups.txt") or die "can not write data into the file RNA_CDS_mpileups.txt: $!";
for my $m (@mfs){
    my ($uuid,$tcga)=split("\\.",$m);
    if ($tcga~~\@cds) {
        print $T $m,"\t",$tcga,"\t","$tcga.bed\n";   
    }   
}
close $T;

#Only run these files absent in the ase_counts dir
`ls $ase_counts > already_done.txt`;

$parsing->vlookup("$ase/RNA_CDS_mpileups.txt",1,"$ase/already_done.txt",1,1,"y","left_ase_grep_NaN.txt");
`grep NaN left_ase_grep_NaN.txt > left_ase_pull_column.txt`;
$parsing->pull_column("left_ase_pull_column.txt","1,2,3","RNA_seq_id_ase_no_tum_norm.txt");

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

opendir (ASE, "$ase/$ase_counts") or die "Can't open directory $ase/$ase_counts: $!\n";
my @asecount = grep{!/^\./ && -f "$ase/$ase_counts/$_"}readdir(ASE);
closedir(ASE);

#Checks if the file contains mpileups that do not have ase_counts
if (-s "RNA_seq_id_ase_no_tum_norm.txt" != 0)
{
            $ase_analysis->compile_ase_no_tum_norm("RNA_seq_id_ase_no_tum_norm.txt","$RNA_Path/$cds_sorted","$mpileups_path",$ase_counts);  
}
#if the above file is empy, it then checks if the ase_counts directory is empty before running the process on all mpileups
elsif(scalar(@asecount) == 0)
{
            $ase_analysis->compile_ase_no_tum_norm("RNA_CDS_mpileups.txt","$RNA_Path/$cds_sorted","$mpileups_path",$ase_counts);
}

#get gene level ase_counts

mkdir "$ase/$gene_level" unless(-d "$ase/$gene_level");
`ls $ase/$gene_level > already_done_genes.txt`;
`ls $mpileups_path/ > mpileups.txt`;

$parsing->vlookup("$ase/mpileups.txt",1,"$ase/already_done_genes.txt",1,1,"y","$ase/left_gene_grep_NaN.txt");
`grep NaN $ase/left_gene_grep_NaN.txt > $ase/left_gene_pull_collumn.txt`;
$parsing->pull_column("$ase/left_gene_pull_collumn.txt",1,"$ase/left_gene_level_ase.txt");

$parsing->vlookup("$ase/RNA_CDS_mpileups.txt",1,"$ase/left_gene_level_ase.txt",1,1,"y","$ase/RNA_seq_id_lookup_grep_non_NaN.txt");
`grep NaN -v $ase/RNA_seq_id_lookup_grep_non_NaN.txt > $ase/RNA_seq_id_lookup_pull_column.txt`;
$parsing->pull_column("$ase/RNA_seq_id_lookup_pull_column.txt","1,2,3","$ase/RNA_seq_id_lookup_comp_gene_faster.txt");

#copy the pbinom.R file
`cp $database_path/pbinom.R $ase`;

#compile_gene_ase_faster will use bed file saved in the directory cds_ase_sorted.
#compile_gene_ase_faster(file that conatins list of ase_counts and associated bed files,cds_sorted_ase directory,path to the refseq.ucsc.ensembl.mrna.hg9.nr.bed file in the Database directory,user defined directory from the command line in script 3.0 or default ase)

if (-s "$ase/RNA_seq_id_lookup_comp_gene_faster.txt" != 0)
{
            $ase_analysis->compile_gene_ase_faster("$ase/RNA_seq_id_lookup_comp_gene_faster.txt","$RNA_Path/$cds_sorted","$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed","$ase","$ase_counts","$gene_level");
}


print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
