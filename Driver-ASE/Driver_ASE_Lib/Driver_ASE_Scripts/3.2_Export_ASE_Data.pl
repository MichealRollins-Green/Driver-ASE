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
    'bsc|b=s' => \my $bad_snps_cnvs, #option to include bad SNPs and CNVs in analysis
    'overlap|o=s' => \my $overlap, #option to overlap RNA-Seq TCGA IDs with matching WGS IDs
    'duplicates|du=s' => \my $duplicates, #option to allow for duplicate ids to seperate any two samples that are associated with one UUID
    'compress|c=s' => \my $compress, #option to compress data
    'remfiles|r=s' => \my $rem_files, #option to remove files when done analysis
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("3.2");

if($help)
{
    $parsing->usage("3.2");
}

if(!defined $disease_abbr)
{
    print "disease was not entered!\n";
    $parsing->usage("3.2");
}

#default to no if no option to includ bad SNPs and CNVs was not specified in the command line
if (!defined $bad_snps_cnvs)
{
    $bad_snps_cnvs = "no";
}
elsif (lc $bad_snps_cnvs ne "y" and lc $bad_snps_cnvs ne "yes" and lc $bad_snps_cnvs ne "n" and lc $bad_snps_cnvs ne "no")
{
    print "Incorrect option specified. Must be yes|y|no|n. \n";
    $parsing->usage("3.2");
}

#default to no if option to overlap ids was not specified in the command line
if (!defined $overlap)
{
    $overlap = "no";
}
elsif (lc $overlap ne "y" and lc $overlap ne "yes" and lc $overlap ne "n" and lc $overlap ne "no")
{
    print "Incorrect option specified. Must be yes|y|no|n. \n";
    $parsing->usage("3.2");
}

#default to yes if option to have duplicate ids is not specified in the command line 
if (!defined $duplicates)
{
    $duplicates = "yes";
}
elsif(lc $duplicates ne "y" and lc $duplicates ne "yes" and lc $duplicates ne "n" and lc $duplicates ne "no")
{
    print "Incorrect option specified. Must be yes|y|no|n. \n";
    $parsing->usage("3.2");
}

#default to no if option to compress files is not specified in the command line
if (!defined $compress)
{
    $compress = "no"
}
elsif(lc $compress ne "y" and lc $compress ne "yes" and lc $compress ne "n" and lc $compress ne "no")
{
    print "Incorrect option specified. Must be yes|y|no|n. \n";
    $parsing->usage("3.2"); 
}

#default to yes if option to remove files is no specified in the command line
if (!defined $rem_files)
{
    $rem_files = "no"
}
elsif(lc $rem_files ne "y" and lc $rem_files ne "yes" and lc $rem_files ne "n" and lc $rem_files ne "no")
{
    print "Incorrect option specified. Must be yes|y|no|n. \n";
    $parsing->usage("3.2");
}

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $tables = "$disease_abbr\_tables";
my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";
my $ase = "$RNA_Path/ase";
my $gene_level = "gene_level";
my $matrix;
my ($bad_cnvs,$bad_snps_bed) = ("bad_cnvs","bad_snps_bed");
my ($problem_snps,$problem_cnvs) = ("problem_snps","problem_cnvs");
my $finished_RNA = "$disease_abbr\_finished_analysis_RNA";
my ($temp,$cds_sorted_ase,$logs) = ("temp","cds_sorted_ase","logs");
my $RNA_table_file = "final_downloadtable_$disease_abbr\_RNA-Seq.txt";
my $WGS_table_file = "final_downloadtable_$disease_abbr\_WGS.txt";
my $Geno_table_file = "$disease_abbr.Genotypes.id2uuid.txt";
my $rna_mpileups = "$ase/rna_mpileups";

#determine what mix of ouput the user wants in their analysis
if((lc $duplicates eq "n" or lc $duplicates eq "n") and (lc $bad_snps_cnvs eq "no" or lc $bad_snps_cnvs eq "n") and (lc $overlap eq "n" or lc $overlap eq "no"))
{
    $matrix = "matrix";
}
elsif((lc $overlap eq "y" or lc $overlap eq "yes") and (lc $duplicates eq "no" or lc $duplicates eq "n") and (lc $bad_snps_cnvs eq "no" or lc $bad_snps_cnvs eq "n"))
{
    $matrix = "matrix_overlap";
}
elsif((lc $overlap eq "y" or lc $overlap eq "yes") and (lc $duplicates eq "yes" or lc $duplicates eq "y") and (lc $bad_snps_cnvs eq "no" or lc $bad_snps_cnvs eq "n"))
{
    $matrix = "matrix_overlap_duplicates";
}
elsif((lc $bad_snps_cnvs eq "yes" or lc $bad_snps_cnvs eq "y") and (lc $overlap eq "n" or lc $overlap eq "no") and (lc $duplicates eq "no" or lc $duplicates eq "n"))
{
    $matrix = "matrix_snps_cnvs";
}
elsif((lc $bad_snps_cnvs eq "yes" or lc $bad_snps_cnvs eq "y") and (lc $duplicates eq "yes" or lc $duplicates eq "y") and (lc $overlap eq "n" or lc $overlap eq "no"))
{
    $matrix = "matrix_snps_cnvs_duplicates";
}
elsif((lc $bad_snps_cnvs eq "yes" or lc $bad_snps_cnvs eq "y") and (lc $overlap eq "y" or lc $overlap eq "yes") and (lc $duplicates eq "no" or lc $duplicates eq "n"))
{
    $matrix = "matrix_snps_cnvs_overlap";
}
elsif((lc $bad_snps_cnvs eq "yes" or lc $bad_snps_cnvs eq "y") and (lc $overlap eq "y" or lc $overlap eq "yes") and (lc $duplicates eq "yes" or lc $duplicates eq "y"))
{
    $matrix = "matrix_snps_cnvs_overlap_duplicates";
}
else
{
    $matrix  = "matrix_with_duplicates";
}

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
    if (!(-e "$Analysispath/$disease_abbr/$tables/$RNA_table_file"))
    {
	print STDERR "$Analysispath/$disease_abbr/$tables/$RNA_table_file does not exist. It was either moved, renamed or deleted,\n";
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

#check if user wants to overlap RNA-Seq,WGS and Genotypes IDs
if (lc $overlap eq "y" || lc $overlap eq "yes")
{
    chdir "$Analysispath/$disease_abbr/$tables/";
    if (!(-f "$Analysispath/$disease_abbr/$tables/$Geno_table_file") and !(-f "$Analysispath/$disease_abbr/$tables/$RNA_table_file") and !(-f "$Analysispath/$disease_abbr/$tables/$WGS_table_file"))
    {
        print "One or more of the table files do not exist: $Geno_table_file, $RNA_table_file, $WGS_table_file\n";
        print "Check the $Analysispath/$disease_abbr/$tables directory and run the table scripts for the tables that may be missing.\n";
        exit;
    }
    else
    {
        open(GENI,"$Geno_table_file") or die "Can't open $Analysispath/$disease_abbr/$tables/$Geno_table_file for input: $!\n";
        open(GN,">$disease_abbr\_Normal_Gen.txt") or die "Can't open $Analysispath/$disease_abbr/$tables/$disease_abbr\_Normal_Gen.txt for output: $!\n";
        open(GT,">$disease_abbr\_Tumor_Gen.txt") or die "Can't open $Analysispath/$disease_abbr/$tables/$disease_abbr\_Tumor_Gen.txt for output: $!\n";
        open(OTH,">$disease_abbr\_Other_Gen.txt") or die "Can't open $Analysispath/$disease_abbr/$tables/$disease_abbr\_Other_Gen.txt for output: $!\n";
        #Format of Genotypes file
        
        #UUID: 1e2b5cf3-a939-4dea-9b77-68faaf871c5a
        #TCGA ID with Sample ID: TCGA-97-7546-01A
        #Sample ID: 1
        #Sample Type: Primary Tumor
        #TCGA ID: TCGA-97-7546
        
        while (my $r = <GENI>)
        {
            chomp($r);
            my @TN = split("\t",$r);
            if ($TN[3] =~ /tumor/ig)
            {
                print GT $r,"\n";
            }
            elsif ($TN[3] =~ /normal/ig)
            {
                print GN $r,"\n";
            }
            else
            {
                print OTH $r,"\n";
            }
        }
        close(GENI);
        close(GN);
        close(GT);
        close OTH;
        
        #Intersect the normal and tumor Genotypes to get matching tumor normal ids
        $ase_analysis->Intersect_Files("$disease_abbr\_Normal_Gen.txt,$disease_abbr\_Tumor_Gen.txt",0,"4,4","$disease_abbr.Genotype_Int.txt");
        
        #parse files to remove non matching tumor/normal
        open(GIN,"$disease_abbr.Genotype_Int.txt") or die "Can't open $Analysispath/$disease_abbr/$tables/$disease_abbr.Genotype_Int.txt for input: $!\n";
        open(GOUT,">$disease_abbr\_Genotype.info.txt") or die "Can't open file $disease_abbr\_Genotype.info.txt: $!\n";
        
        while (my $r = <GIN>)
        {
            chomp($r);
            my @Genos = split("\t",$r);
            
            if ($Genos[1] =~ /normal/ig && $Genos[1] =~ /tumor/ig)
            {
                print GOUT $r,"\n";
            }
        }
        close(GIN);
        close(GOUT);
        
        #Intersect parsed file with RNA-Seq and WGS
        $ase_analysis->Intersect_Files("$RNA_table_file,$WGS_table_file,$disease_abbr\_Genotype.info.txt",0,"1,1,2","$disease_abbr\_Intersect.txt");
        
        #parse intersected file to remove IDs that do not intersect with RNA-Seq, WGS and Genotypes
        open(INT,"$disease_abbr\_Intersect.txt") or die "can't open $Analysispath/$disease_abbr/$tables/$disease_abbr\_Intersect.txt for input: $!\n";
        open(INTOUT,">$disease_abbr\_RNA_WGS_Geno.txt") or die "Can't open file $Analysispath/$disease_abbr/$tables/$disease_abbr\_RNA_WGS_Geno.txt: $!\n";
        
        while (my $r = <INT>)
        {
            chomp($r);
            my @RWG = split("\t",$r);
            
            if ($RWG[1] =~ /RNA/ig && $RWG[1] =~ /WGS/ig && $RWG[1] =~ /Genotype/ig)
            {
                print INTOUT $r,"\n";
            }
        }
        close(INT);
        close(INTOUT);  
    }
    chdir $RNA_Path;
}

 #get gene level ase files
print "Opening directory $ase/$gene_level for reading.\n";
opendir (GENE,"$ase/$gene_level");
open (my $O,">$ase/gene_level_ases.txt") or die "can not write data into the file $ase/gene_level_ases.txt: $!";

my @ases = grep{!/^\./ && -f "$ase/$gene_level/$_"} readdir(GENE);
closedir(GENE);
@ases = map{
    my $l = $_;
    chomp($l);
    my @as = split("\\.",$l);
    my $id = $as[-1];
    $id =~ s/-[0-9]+[a-zA-Z]$//;
    print $O join("\t",@as),"\t$l\t$id\n";
}@ases;
close($O);

#parse files to get the ones with overlappng RNA-Seq, WGS and Genotypes IDs
if (lc $overlap eq "y" || lc $overlap eq "yes")
{
    $parsing->vlookup("$ase/gene_level_ases.txt",4,"$Analysispath/$disease_abbr/$tables/$disease_abbr\_RNA_WGS_Geno.txt",3,3,"y","$ase/gene_NaN.txt");
    `grep -v -i NaN $ase/gene_NaN.txt > $ase/gene_pull.txt`;
    $parsing->pull_column("$ase/gene_pull.txt","1,2,3","$ase/gene_level_ases.txt");
}

#check if user is using bad SNPs and CNVs in analysis
if (lc $bad_snps_cnvs eq "y" || lc $bad_snps_cnvs eq "yes")
{           
    $parsing->vlookup("$ase/gene_level_ases.txt",3,"$ase/RNA_CDS_mpileups.txt",1,2,"y","$ase/All_ASEs_grep_non_NaN.txt");
    `grep NaN -v -i $ase/All_ASEs_grep_non_NaN.txt > $ase/All_ASEs.txt`;
    
    $parsing->vlookup("$ase/All_ASEs.txt",1,"$Analysispath/$disease_abbr/$tables/final_downloadtable_$disease_abbr\_RNA-Seq.txt",1,"6,3","y","$ase/ASEs4lookup_pull_column.txt");
    
    #if overlap is not specified, extra columns will be present in the file making the columns to pull different.
    if (lc $overlap eq "n" or lc $overlap eq "no")
    {
        $parsing->pull_column("$ase/ASEs4lookup_pull_column.txt","3,2,6,7","$ase/ASEs4lookup_sort.txt");
    }
    else
    {
        $parsing->pull_column("$ase/ASEs4lookup_pull_column.txt","3,2,5,6","$ase/ASEs4lookup_sort.txt");
    }
    
    `sort -k 2,2 -k 4,4n $ase/ASEs4lookup_sort.txt > $ase/ASEs4lookup.txt`;
    
    #pull_TN acts as a SAS lag function or retain function
    #pull_TN(ase bad cnv file, bad_cnvs directory,bad_snps_bed directory,output file)
    $ase_analysis->pull_TN("$ase/ASEs4lookup.txt","$RNA_Path/$bad_cnvs","$RNA_Path/$bad_snps_bed","$RNA_Path/have_cnv_snps.txt",$disease_abbr);
    
    #functionality of have_cnv_snp is to keep IDs having both bad_cnv and bad_snps
    #have_cnv_snp(output file created in pull_TN,output file)
    $ase_analysis->have_cnv_snp("$RNA_Path/have_cnv_snps.txt","$RNA_Path/have_cnv_snps_sort.txt");
    
    `sort -k 2,2 $RNA_Path/have_cnv_snps_sort.txt > $RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt`;

    #Now link Tumor and Normals gene_level files with TCGA IDs (12 chars)
    #Output files will be tumor_ff and normal_ff
    #mk_files(the TN snp cnv lookup file created from the sort command, path to the user defined directory from the command line in script 3.0 or default ase)
    $ase_analysis->mk_files_snps_cnvs("$RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt","$ase","$gene_level",$duplicates);
}
else
{
    #if overlap is no then there will be an additional column in the file that is not there when overlap is yes. This column is removed
    if (lc $overlap eq "n" or lc $overlap eq "no")
    {
        $parsing->vlookup("$ase/gene_level_ases.txt",1,"$Analysispath/$disease_abbr/$tables/$RNA_table_file",1,"2,3,4,5,6,7","y","$ase/ASEs4lookup_pull.txt");
        $parsing->pull_column("$ase/ASEs4lookup_pull.txt","1,2,3,4,6,7,8,9,10","$ase/ASEs4lookup.txt") unless (lc $overlap eq "y" or lc $overlap eq "yes");
    }
    else
    {
        $parsing->vlookup("$ase/gene_level_ases.txt",1,"$Analysispath/$disease_abbr/$tables/$RNA_table_file",1,"2,3,4,5,6,7","y","$ase/ASEs4lookup.txt");
    }
    #Output files will be tumor_ff and normal_ff
    $ase_analysis->mk_files_tum_norm("$ase/ASEs4lookup.txt","$ase/$gene_level",$duplicates);
}

mkdir "$ase/$matrix" unless(-d "$matrix");

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

#check if user wants duplicated ids
if (lc $duplicates eq "y" || lc $duplicates eq "yes")
{
    `sort -k 1,1 collabels.txt > collabels_sort.txt`;
}
else
{
    `sort -u -k 1,1 collabels.txt > collabels_sort.txt`;
}

`mv $RNA_Path/collabels_sort.txt $ase/$matrix/tumor.collabels`;

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

#check if user wants duplicated ids
if (lc $duplicates eq "y" || lc $duplicates eq "yes")
{
    `sort -k 1,1 collabels.txt > collabels_sort.txt`;
}
else
{
    `sort -u -k 1,1 collabels.txt > collabels_sort.txt`;
}

`mv $RNA_Path/collabels_sort.txt $ase/$matrix/normal.collabels`;
 
$parsing->vlookup("$ase/$matrix/rowlabels.txt",1,"$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",4,"1,2,3,4,5,6,7,8,9,10,11,12","n","$ase/$matrix/rowlabels.bed");

if (lc $bad_snps_cnvs eq "y" || lc $bad_snps_cnvs eq "yes")
{
    #problem zones -> pull bed from rowlabels and make lists of genes
    #rowlabels.bed will be used by script run_problem_cnvs and run_problem_snps!
    #run_problem_cnvs(bad_cnvs directory,the TN snp cnv lookup file created from the sort command,path to RNA_Seq_Analysis directory,matrix directory)
    $ase_analysis->run_problem_cnvs("$RNA_Path/$bad_cnvs","$RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt","$ase","$matrix","$problem_snps");
    #run_problem_snps(bad_snps_bed directory,the TN snp cnv lookup file created from the sort command,path to RNA_Seq_Analysis directory,matrix directory)
    $ase_analysis->run_problem_snps("$RNA_Path/$bad_snps_bed","$RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt","$ase","$matrix","$problem_cnvs");
    
    open(TUM,"tumor_ff") or die "Can't open tumor_ff: $!\n";
    open(CNV,">cnv_ff") or die "Can't open cnv_ff: $!\n";
    
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
}

#check if user wants files used in analysis deleted
if (lc $rem_files eq "y" or lc $rem_files eq "yes")
{
    my @del_files = $parsing->get_only_files_in_dir("$RNA_Path");

    @del_files = grep{!/\.bim/ and !/\.fam/ and !/\.log/ and !/\.bed/}@del_files;
    #Delete files in the RNA_Seq_Analysis directory
    for(my $i = 0;$i < scalar(@del_files);$i++)
    {
	`rm "$del_files[$i]"`;
    }
    
    chdir $ase;
    
    undef @del_files;
    @del_files = $parsing->get_only_files_in_dir("$ase");
    #Delete files in the ase directory
    for(my $i = 0;$i < scalar(@del_files);$i++)
    {
	`rm "$del_files[$i]"`;
    }
}

chdir $RNA_Path;

#check if user wants data to be compressed
if (lc $compress eq "y" or lc $compress eq "yes")
{
    if (-e "$disease_abbr\_ASE_Rst.tar.gz")
    {
	#check if user included bad SNPs and CNVs in analysis
        if (-d $bad_snps_bed and -d $bad_cnvs)
        {
            print "Now running tar -zcvkf $disease_abbr\_ASE_Rst.tar.gz --exclude \"$ase/$temp\" --exclude \"$ase/$logs\" $bad_snps_bed $bad_cnvs $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*.\n";
            
	    `tar -zcvkf $disease_abbr\_ASE_Rst.tar.gz --exclude "$ase/$temp" --exclude "$ase/$logs" $bad_snps_bed $bad_cnvs $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*`;
        }
        else
        {
	    print "Now running tar -zcvkf $disease_abbr\_ASE_Rst.tar.gz --exclude \"$ase/$temp\" --exclude \"$ase/$logs\" $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*.\n";
            
            `tar -zcvkf $disease_abbr\_ASE_Rst.tar.gz --exclude "$ase/$temp" --exclude "$ase/$logs" $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*`;
        }
    }
    else
    {
	#check if user included bad SNPs and CNVs in analysis
        if (-d $bad_snps_bed and -d $bad_cnvs)
        {
        print "Now running tar -zcvf $disease_abbr\_ASE_Rst.tar.gz --exclude \"$ase/$temp\" --exclude \"$ase/$logs\" $bad_snps_bed $bad_cnvs $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*.\n";
            
        `tar -zcvf $disease_abbr\_ASE_Rst.tar.gz --exclude "$ase/$temp" --exclude "$ase/$logs" $bad_snps_bed $bad_cnvs $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*`;
        }
	else
        {
            print "Now running tar -zcvf $disease_abbr\_ASE_Rst.tar.gz --exclude \"$ase/$temp\" --exclude \"$ase/$logs\" $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*.\n";
            
            `tar -zcvf $disease_abbr\_ASE_Rst.tar.gz --exclude "$ase/$temp" --exclude "$ase/$logs" $ase/* $disease_abbr\_TN_TCGA_Imputation.* $disease_abbr\_TN_TCGA_All.*`;
        }
    }

    mkdir "$Analysispath/$disease_abbr/$finished_RNA" unless(-d "$Analysispath/$disease_abbr/$finished_RNA");

    `mv $disease_abbr\_ASE_Rst.tar.gz $Analysispath/$disease_abbr/$finished_RNA`;
}

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
