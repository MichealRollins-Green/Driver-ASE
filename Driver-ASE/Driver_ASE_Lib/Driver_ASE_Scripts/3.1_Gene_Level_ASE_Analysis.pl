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
use autodie;
no warnings 'once';

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $ase_analysis = Driver_ASE_Lib::Gene_Level_ASE_Analysis->new;
my $parsing = Driver_ASE_Lib::Parsing_Routines->new;

GetOptions(
    'cancer|c=s' => \my $cancer_type,#e.g. OV
    'overlap|o=s' => \my $overlap, #option to overlap RNA-Seq TCGA IDs with matching WGS IDs
    'overlapSelect|S=s' => \my $overlapSelect,#path to overlapSelect
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("3.1");

if ($help)
{
    $parsing->usage("3.1");
}

if (!defined $cancer_type)
{
    print STDERR "Cancer type was not entered!\n";
    $parsing->usage("3.1");
}

#default to no if option to overlap ids was not specified in the command line
if (!defined $overlap)
{
    $overlap = "no";
}
elsif (lc $overlap ne "y" and lc $overlap ne "yes" and lc $overlap ne "n" and lc $overlap ne "no")
{
    print STDERR "Incorrect option specified! Must be yes|y|no|n.\n";
    $parsing->usage("3.1");
}

if (!defined $overlapSelect)
{
    $overlapSelect = "overlapSelect";
}

$parsing->CheckSoftware("$overlapSelect");

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $tables = "$cancer_type\_tables";
my $RNA_Path = "$Analysispath/$cancer_type/RNA_Seq_Analysis";
my $ase = "$RNA_Path/ase";
my $mpileups_dir = "rna_mpileups";
my ($logs,$ase_counts,$cds_sorted,$gene_level) = ("logs","ase_counts","cds_sorted","gene_level");
my $refseq = "refseq.ucsc.ensembl.mrna.hg9.nr.bed";
my ($RNA_CDS_Mpileups,$RNA_CDS_Mpileups_overlap) = ("RNA_CDS_mpileups.txt","RNA_CDS_mpileups_overlap.txt");
my ($RNA_table_file,$WGS_table_file,$Geno_table_file) = ("final_downloadtable_$cancer_type\_RNA-Seq.txt","final_downloadtable_$cancer_type\_WGS.txt","$cancer_type.Genotypes.id2uuid.txt");
my ($RNA_table_overlap,$WGS_table_overlap,$Geno_table_overlap) = ("final_downloadtable_$cancer_type\_RNA-Seq_overlap.txt","final_downloadtable_$cancer_type\_WGS_overlap.txt","$cancer_type.Genotypes.id2uuid_overlap.txt");
my $Intersect = "$cancer_type\_RNA_WGS_Geno.txt";
my $rna_mp_file = "rna_mpileups_overlap.txt";
my $Exp_Strategy = "RNA-Seq";
my $finished_RNA = "$cancer_type\_finished_analysis_RNA";
my ($ase_counts_done,$gene_level_done) = ("$ase/ase_counts_done.txt","$ase/gene_level_done.txt");

my @mfs;

$parsing->check_directory_existence("$database_path","$Analysispath","$Analysispath/$cancer_type","$RNA_Path","$RNA_Path/$cds_sorted","$ase","$ase/$mpileups_dir"); #check if directories or files exist

$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

chdir "$ase";

mkdir "$logs" unless(-d "$logs");
mkdir "$ase_counts" unless(-d "$ase_counts");

my $cdsmp;
#check if user wants to overlap RNA-Seq,WGS and Genotypes IDs
if (lc $overlap eq "y" || lc $overlap eq "yes")
{
    $parsing->check_directory_existence("$Analysispath/$cancer_type/$tables/$Geno_table_file","$Analysispath/$cancer_type/$tables/$RNA_table_file","$Analysispath/$cancer_type/$tables/$WGS_table_file");
    
    if (!(-s "$Analysispath/$cancer_type/$tables/$RNA_table_file" == 0) and !(-s "$Analysispath/$cancer_type/$tables/$WGS_table_file" == 0) and !(-s "$Analysispath/$cancer_type/$tables/$Geno_table_file" == 0))
    {
        $parsing->Overlap_RNA_WGS_Geno("$Analysispath/$cancer_type/$tables","$RNA_table_file","$WGS_table_file","$Geno_table_file","$RNA_table_overlap","$WGS_table_overlap","$Geno_table_overlap","$Intersect","$Exp_Strategy","$cancer_type");
       
        chdir "$ase"; #changes back to this directory as Overlap_RNA_WGS_Geno changes to the table directory for the cancer type
     
        my @rmp = `ls $ase/$mpileups_dir`;
        
        #fromat in files must be UUID.TCGA_ID
        $ase_analysis->ase_filter_overlap(\@rmp,"$Analysispath/$cancer_type/$tables/$RNA_table_overlap",$rna_mp_file);
        
        @mfs = `cat $rna_mp_file`;
        
        #match cds with its corresponding mpileup file;
        open ($cdsmp, ">$ase/$RNA_CDS_Mpileups_overlap");
    }
}
else
{
    #Get mpileup files;
    opendir (MP, "$ase/$mpileups_dir");
    @mfs = grep{!/^\./ && -f "$ase/$mpileups_dir/$_"}readdir(MP);
    closedir (MP);
    
    #match cds with its corresponding mpileup file;
    open ($cdsmp, ">$ase/$RNA_CDS_Mpileups");
}

#Get cds files;
opendir (CD, "$RNA_Path/$cds_sorted");
my @cds= grep{!/^\./ && -f "$RNA_Path/$cds_sorted/$_"} readdir(CD);
closedir (CD);
@cds = map{
s/\.bed//;$_;
}@cds;

for my $m (@mfs)
{
    chomp($m);
    my ($uuid,$tcga) = split("\\.",$m);
    if ($tcga~~\@cds)
    {
        print $cdsmp $m,"\t",$tcga,"\t","$tcga.bed\n";
    }
}
close ($cdsmp);

#Only run these files absent in the $ase_counts directory
if (lc $overlap eq "y" or lc $overlap eq "yes")
{
    my @done_ase = `ls $ase_counts`;
    
    #fromat in files must be UUID.TCGA_ID
    $ase_analysis->ase_filter_overlap(\@done_ase,"$Analysispath/$cancer_type/$tables/$RNA_table_overlap",$ase_counts_done);
    $parsing->vlookup("$ase/$RNA_CDS_Mpileups_overlap",1,"$ase_counts_done",1,1,"y","left_ase_grep_NaN.txt");
}
else
{
    `ls $ase_counts > $ase_counts_done`;
    $parsing->vlookup("$ase/$RNA_CDS_Mpileups",1,"$ase_counts_done",1,1,"y","left_ase_grep_NaN.txt");
}

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

opendir (ASE, "$ase/$ase_counts");
my @asecount = grep{!/^\./ && -f "$ase/$ase_counts/$_"}readdir(ASE);
closedir (ASE);

#compile_ase_no_tum_norm(file that conatins list of mpileups and associated bed files, $RNA_Path/$cds_sorted (path to directory where cds sorted beds are stored), $mpileups_dir (path to the directory where RNA-Seq mpileups are strored), $ase_counts (path the the directory where ASE count data will be stored))
#Checks if the file contains mpileups that do not have ase_counts
if (-s "RNA_seq_id_ase_no_tum_norm.txt" != 0)
{
    $ase_analysis->compile_ase_no_tum_norm("RNA_seq_id_ase_no_tum_norm.txt","$RNA_Path/$cds_sorted","$ase/$mpileups_dir",$ase_counts);
}
#if the above file is empy, it then checks if the ase_counts directory is empty before running the process on all mpileups
elsif (scalar(@asecount) == 0)
{
    if (lc $overlap eq "y" or lc $overlap eq "yes")
    {
        $ase_analysis->compile_ase_no_tum_norm("$RNA_CDS_Mpileups_overlap","$RNA_Path/$cds_sorted","$ase/$mpileups_dir",$ase_counts);
    }
    else
    {
        $ase_analysis->compile_ase_no_tum_norm("$RNA_CDS_Mpileups","$RNA_Path/$cds_sorted","$ase/$mpileups_dir",$ase_counts);
    }
}

mkdir "$ase/$gene_level" unless(-d "$ase/$gene_level");

#Only run these files absent in the $ase_counts directory
if (lc $overlap eq "y" or lc $overlap eq "yes")
{
    my @done_gene = `ls $ase/$gene_level`;
 
    #fromat in files must be UUID.TCGA_ID
    $ase_analysis->ase_filter_overlap(\@done_gene,"$Analysispath/$cancer_type/$tables/$RNA_table_overlap","$gene_level_done");
    #$parsing->vlookup("$ase/$RNA_CDS_Mpileups_overlap",1,"$gene_level_done",1,1,"y","left_ase_grep_NaN.txt");
    $parsing->vlookup("$ase/$RNA_CDS_Mpileups_overlap",1,"$gene_level_done",1,1,"y","$ase/left_gene_grep_NaN.txt");
}
else
{
    `ls $ase/$gene_level > $gene_level_done`;
    #$parsing->vlookup("$ase/mpileups.txt",1,"$gene_level_done",1,1,"y","$ase/left_gene_grep_NaN.txt");
    $parsing->vlookup("$ase/$RNA_CDS_Mpileups",1,"$gene_level_done",1,1,"y","left_gene_grep_NaN.txt");
}

`grep NaN $ase/left_gene_grep_NaN.txt > $ase/left_gene_pull_collumn.txt`;
$parsing->pull_column("$ase/left_gene_pull_collumn.txt",1,"$ase/left_gene_level_ase.txt");

if (lc $overlap eq "y" or lc $overlap eq "yes")
{
    $parsing->vlookup("$ase/$RNA_CDS_Mpileups_overlap",1,"$ase/left_gene_level_ase.txt",1,1,"y","$ase/RNA_seq_id_lookup_grep_non_NaN.txt");
}
else
{
    $parsing->vlookup("$ase/$RNA_CDS_Mpileups",1,"$ase/left_gene_level_ase.txt",1,1,"y","$ase/RNA_seq_id_lookup_grep_non_NaN.txt");
}

`grep NaN -v $ase/RNA_seq_id_lookup_grep_non_NaN.txt > $ase/RNA_seq_id_lookup_pull_column.txt`;
$parsing->pull_column("$ase/RNA_seq_id_lookup_pull_column.txt","1,2,3","$ase/RNA_seq_id_lookup_comp_gene_faster.txt");

#copy the pbinom.R file
`cp $database_path/pbinom.R $ase`;

#compile_gene_ase_faster($ase/RNA_seq_id_lookup_comp_gene_faster.txt (file that conatins list of ASE count files and associated bed files), $RNA_Path/$cds_sorted (path the the directory that contains cds sorted files), $database_path/$refseq (path to the $refseq file in the $database_path directory), $ase (path to directory where ASE analysis  data is stored), $ase_counts (directory where ase count data was stored), $gene_level (directory where gene level ase data will be stored), $overlapSelect (overlapSelect command))

if (-s "$ase/RNA_seq_id_lookup_comp_gene_faster.txt" != 0)
{
    $ase_analysis->compile_gene_ase_faster("$ase/RNA_seq_id_lookup_comp_gene_faster.txt","$RNA_Path/$cds_sorted","$database_path/$refseq","$ase","$ase_counts","$gene_level","$overlapSelect");
}

print "All jobs have finished for $cancer_type.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
