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
    'overlapSelect|S=s' => \my $overlapSelect, #path to overlapSelect
    'badsnpcnvs|b=s' => \my $bad_snps_cnvs, #option to include bad SNPs and CNVs in analysis
    'overlap|o=s' => \my $overlap, #overlap option. only enter in if option was y|yes in script 3.1
    'duplicates|d=s' => \my $duplicates, #option to allow for duplicate ids to seperate any two samples that are associated with one UUID
    'remfiles|r=s' => \my $rem_files, #option to remove files when done analysis
    'archive|a=s' => \my $archive, #option to archive files
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("3.2");

if($help)
{
    $parsing->usage("3.2");
}

if(!defined $cancer_type)
{
    print STDERR "Cancer type was not entered!\n";
    $parsing->usage("3.2");
}

#default to no if no option to includ bad SNPs and CNVs was not specified in the command line
if (!defined $bad_snps_cnvs)
{
    $bad_snps_cnvs = "no";
}
elsif (lc $bad_snps_cnvs ne "y" and lc $bad_snps_cnvs ne "yes" and lc $bad_snps_cnvs ne "n" and lc $bad_snps_cnvs ne "no")
{
    print STDERR "Incorrect option specified! Must be yes|y|no|n.\n";
    $parsing->usage("3.2");
}
else
{
    if (!defined $overlapSelect)
    {
        $overlapSelect = "overlapSelect";
    }
    $parsing->CheckSoftware("$overlapSelect");
}

#default to no if option to overlap ids was not specified in the command line
if (!defined $overlap)
{
    $overlap = "no";
}
elsif (lc $overlap ne "y" and lc $overlap ne "yes" and lc $overlap ne "n" and lc $overlap ne "no")
{
    print STDERR "Incorrect option specified! Must be yes|y|no|n.\n";
    $parsing->usage("3.2");
}

#default to yes if option to have duplicate ids is not specified in the command line 
if (!defined $duplicates)
{
    $duplicates = "yes";
}
elsif(lc $duplicates ne "y" and lc $duplicates ne "yes" and lc $duplicates ne "n" and lc $duplicates ne "no")
{
    print STDERR "Incorrect option specified! Must be yes|y|no|n.\n";
    $parsing->usage("3.2");
}

#default to yes if option to remove files is no specified in the command line
if (!defined $rem_files)
{
    $rem_files = "no"
}
elsif(lc $rem_files ne "y" and lc $rem_files ne "yes" and lc $rem_files ne "n" and lc $rem_files ne "no")
{
    print STDERR "Incorrect option specified! Must be yes|y|no|n.\n";
    $parsing->usage("3.2");
}

if (!defined $archive)
{
    $archive = "no";
}
elsif(lc $archive ne "y" and lc $archive ne "yes" and lc $archive ne "n" and lc $archive ne "no")
{
    print STDERR "Incorrect option specified! Must be yes|y|no|n.\n";
    $parsing->usage("3.2");    
}

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $tables = "$cancer_type\_tables";
my $RNA_Path = "$Analysispath/$cancer_type/RNA_Seq_Analysis";
my $ase = "$RNA_Path/ase";
my ($maps,$peds,$phased,$cds_sorted,$cds_plink,$ase_counts,$gene_level,$matrix,$mpileups_dir) = ("maps","peds","phased","cds_sorted","cds_plink","ase_counts","gene_level","matrix","rna_mpileups");
my ($bad_cnvs,$bad_snps_bed) = ("bad_cnvs","bad_snps_bed");
my ($problem_snps,$problem_cnvs) = ("problem_snps","problem_cnvs");
my $finished_RNA = "$cancer_type\_finished_analysis_RNA";
my ($temp,$logs) = ("temp","logs");
my ($RNA_table_file,$RNA_table_overlap) = ("final_downloadtable_$cancer_type\_RNA-Seq.txt","final_downloadtable_$cancer_type\_RNA-Seq_overlap.txt");
my ($RNA_CDS_Mpileups,$RNA_CDS_Mpileups_overlap) = ("RNA_CDS_mpileups.txt","RNA_CDS_mpileups_overlap.txt");
my $gene_level_ases = "gene_level_ases.txt";
my ($ase_compress_file,$impute_compress_file) = ("$cancer_type\_ASE_Rst","$cancer_type\_Imputation_Haps.tar.gz");
my ($ase_counts_archive,$gene_level_archive,$cds_sorted_archive,$RNA_CDS_Mpileups_archive) = ("$ase_counts\_to_archive.txt","$gene_level\_to_archive.txt","$cds_sorted\_to_archive.txt","RNA_CDS_mpileups_to_archive.txt");
my ($maps_archive,$peds_archive,$phased_archive,$bad_snp_archive,$bad_cnv_archive,$problem_snp_archive,$problem_cnv_archive,$Impute_archive,$matrix_archive,$impute_haps_archive) = ("$maps\_to_archive.txt","$peds\_to_archive.txt","$phased\_to_archive.txt","$bad_snps_bed\_to_archive.txt","$bad_cnvs\_to_archive.txt","$problem_snps\_to_archive.txt","$problem_cnvs\_to_archive.txt","Imputation_to_archive.txt","$matrix\_to_archive.txt","$cds_plink\_Imputation_Haps_to_archive.txt");

$parsing->check_directory_existence("$database_path","$Analysispath","$Analysispath/$cancer_type","$Analysispath/$cancer_type/$tables","$Analysispath/$cancer_type/$tables/$RNA_table_file","$RNA_Path","$ase","$ase/$gene_level"); #check if directories or files exist
$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

 chdir $ase;

my @ases;
if (lc $overlap eq "y"  or lc $overlap eq "yes")
{
    #get gene level ase files
    print "Getting matching overlapped files $ase/$gene_level.\n";
    `ls $ase/$gene_level > gene_level_look.txt`;
    $parsing->vlookup("gene_level_look.txt",1,"$RNA_CDS_Mpileups_overlap",1,1,"y","gene_level_NaN.txt");
    my @gen_nan = `cat gene_level_NaN.txt`;
    open(GENO,">gene_level_overlap.txt");
    
    for my $gen (@gen_nan)
    {
        chomp($gen);
        if ($gen !~ /\tNaN+$/i)
        {
            my $split_gen = [split("\t",$gen)]->[0]; #There are 2 columns but only the first one is needed
            print GENO $split_gen,"\n";
        }            
    }
    close (GENO);
    @ases = `cat gene_level_overlap.txt`;
}
else
{
    #get gene level ase files
    print "Opening directory $ase/$gene_level for reading.\n";
    opendir (GENE,"$ase/$gene_level");
    @ases = grep{!/^\./ && -f "$ase/$gene_level/$_"} readdir(GENE);
    closedir (GENE);           
}

open (my $O,">$ase/$gene_level_ases");

@ases = map{
    my $l = $_;
    chomp($l);
    my @as = split("\\.",$l);
    my $id = $as[-1];
    $id =~ s/-[0-9]+[a-zA-Z]$//;
    print $O join("\t",@as),"\t$l\t$id\n";
}@ases;
close ($O);

#check if user is using bad SNPs and CNVs in analysis
if (lc $bad_snps_cnvs eq "y" || lc $bad_snps_cnvs eq "yes")
{
    $parsing->check_directory_existence("$RNA_Path/$bad_cnvs","$RNA_Path/$bad_snps_bed");
    
    if (lc $overlap eq "y" or lc $overlap eq "yes")
    {
        $parsing->vlookup("$ase/$gene_level_ases",3,"$ase/$RNA_CDS_Mpileups_overlap",1,2,"y","$ase/All_ASEs_grep_non_NaN.txt");
    }
    else
    {
        $parsing->vlookup("$ase/$gene_level_ases",3,"$ase/$RNA_CDS_Mpileups",1,2,"y","$ase/All_ASEs_grep_non_NaN.txt");
    }
  
    `grep -v NaN $ase/All_ASEs_grep_non_NaN.txt > $ase/All_ASEs.txt`;
    
    $parsing->vlookup("$ase/All_ASEs.txt",1,"$Analysispath/$cancer_type/$tables/$RNA_table_file",1,"6,3","y","$ase/ASEs4lookup_pull_column.txt");
    
    $parsing->pull_column("$ase/ASEs4lookup_pull_column.txt","3,2,6,7","$ase/ASEs4lookup_sort.txt");

    `sort -k 2,2 -k 4,4n $ase/ASEs4lookup_sort.txt > $ase/ASEs4lookup.txt`;

    #pull_TN acts as a SAS lag function or retain function
    #pull_TN($ase/ASEs4lookup.txt file that contains bad cnvs/snps, $RNA_Path/$bad_cnvs (path to directory that holds bad cnv files), $RNA_Path/$bad_snps_bed (path to directory that holds bad snp bed files),output file, $cancer_type (cancer type e.g. PRAD))
    $ase_analysis->pull_TN("$ase/ASEs4lookup.txt","$RNA_Path/$bad_cnvs","$RNA_Path/$bad_snps_bed","$RNA_Path/have_cnv_snps.txt",$cancer_type);

    #functionality of have_cnv_snp is to keep IDs having both bad_cnv and bad_snps
    #have_cnv_snp($RNA_Path/have_cnv_snps.txt (output file created in pull_TN), $RNA_Path/have_cnv_snps_sort.txt (output file))
    $ase_analysis->have_cnv_snp("$RNA_Path/have_cnv_snps.txt","$RNA_Path/have_cnv_snps_sort.txt");
    
    `sort -k 2,2 $RNA_Path/have_cnv_snps_sort.txt > $RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt`;

    #Now link Tumor and Normals gene_level files with TCGA IDs (12 chars)
    #Output files will be tumor_ff and normal_ff
    #mk_files($RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt (the TN snp cnv lookup file created from the sort command), $ase (path to directory where ASE analysis data is stored), $gene_level (directory where gene level files are stored), $duplicates (option for duplicates))
    $ase_analysis->mk_files_snps_cnvs("$RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt","$ase","$gene_level",$duplicates);
}
else
{
    #if overlap is no then there will be an additional column in the file that is not there when overlap is yes. This column is removed
    if (lc $overlap eq "n" or lc $overlap eq "no")
    {
        $parsing->vlookup("$ase/$gene_level_ases",1,"$Analysispath/$cancer_type/$tables/$RNA_table_overlap",1,"3,4,5,6,7","y","$ase/ASEs4lookup.txt");
        #$parsing->pull_column("$ase/ASEs4lookup_pull.txt","1,2,3,4,6,7,8,9,10","$ase/ASEs4lookup.txt") unless (lc $overlap eq "y" or lc $overlap eq "yes");
    }
    else
    {
        $parsing->vlookup("$ase/$gene_level_ases",1,"$Analysispath/$cancer_type/$tables/$RNA_table_file",1,"3,4,5,6,7","y","$ase/ASEs4lookup.txt");
    }

    #Output files will be tumor_ff and normal_ff
    $ase_analysis->mk_files_tum_norm("$ase/ASEs4lookup.txt","$ase/$gene_level",$duplicates);
}

mkdir "$ase/$matrix" unless(-d "$ase/$matrix");
`rm -rf $ase/$matrix/*`;

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
    `sort -k 1,1 $RNA_Path/collabels.txt > collabels_sort.txt`;
}
else
{
    `sort -u -k 1,1 $RNA_Path/collabels.txt > collabels_sort.txt`;
}

`mv collabels_sort.txt $ase/$matrix/tumor.collabels`;

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
    `sort -k 1,1 $RNA_Path/collabels.txt > collabels_sort.txt`;
}
else
{
    `sort -u -k 1,1 $RNA_Path/collabels.txt > collabels_sort.txt`;
}

`mv collabels_sort.txt $ase/$matrix/normal.collabels`;
 
$parsing->vlookup("$ase/$matrix/rowlabels.txt",1,"$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",4,"1,2,3,4,5,6,7,8,9,10,11,12","n","$ase/$matrix/rowlabels.bed");

if (lc $bad_snps_cnvs eq "y" || lc $bad_snps_cnvs eq "yes")
{
    #problem zones -> pull bed from rowlabels and make lists of genes
    #rowlabels.bed will be used by script run_problem_cnvs and run_problem_snps!
    #run_problem_cnvs($RNA_Path/$bad_cnvs (directory that holds bad cnv files), $RNA_Path/RNA_seq_bad_snp_cnv_TN_look.tx (the TN snp cnv lookup file created from the sort command), $ase (path to directory that holds ase analysis results), $matrix (directory that holds matrix data), $problem_cnvs (directory that will store problem cnv data), $overlapSelect (overlapSelect command))
    $ase_analysis->run_problem_cnvs("$RNA_Path/$bad_cnvs","$RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt","$ase","$matrix","$problem_cnvs","$overlapSelect");
    #run_problem_snps($RNA_Path/$bad_snps_bed (directory that holds bad snp bed files), $RNA_Path/RNA_seq_bad_snp_cnv_TN_look.tx (the TN snp cnv lookup file created from the sort command), $ase (path to directory that holds ase analysis results), $matrix (directory that holds matrix data), $problem_snps (directory that will store problem snp data), $overlapSelect (overlapSelect command))
    $ase_analysis->run_problem_snps("$RNA_Path/$bad_snps_bed","$RNA_Path/RNA_seq_bad_snp_cnv_TN_look.txt","$ase","$matrix","$problem_snps","$overlapSelect");
    
   open (TUM,"tumor_ff");
   open (CNV,">cnv_ff");
    
    #changes gene_level in path name to problem_cnvs
    while (my $r = <TUM>)
    {
        chomp($r);
        my @a = split("\t",$r);
        my $pcnv_path = $a[0];
        $pcnv_path =~ s/$gene_level//;
        $a[0] =~ s/$gene_level/$problem_cnvs/;
        print CNV join("\t",@a), "\n";
    }
    close (TUM);
    close (CNV);

    #matricize_version_two(index_file,print_file,index_column,print_index_column,print_output_column)
    $parsing->matricize_version_two("tumor_ff","cnv_ff",4,1,2);
    
    `mv matrix.tab $ase/$matrix/problem.cnv`;
    
   open (TUM,"tumor_ff");
   open (SNP,">snp_ff");
    
    #changes gene_level in path name to problem_snps
    while (my $r = <TUM>)
    {
        chomp($r);
        my @a = split("\t",$r);
        my $psnp_path = $a[0];
        $psnp_path =~ s/$gene_level//;
        $a[0] =~ s/$gene_level/$problem_snps/;
        print SNP join("\t",@a), "\n";
    }
    close (SNP);
    close (TUM);
    
    $parsing->matricize_version_two("tumor_ff","snp_ff",4,1,2);
    
    `mv matrix.tab $ase/$matrix/problem.snp`;
}

#check if user wants files used in analysis deleted
if (lc $rem_files eq "y" or lc $rem_files eq "yes")
{
    chdir $RNA_Path;
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

if (lc $archive eq "y" or lc $archive eq "yes")
{
    #determine what the name of the archive will be based on what data the user include in their analysis
    if (lc $duplicates eq "y" or lc $duplicates eq "yes")
    {
        $ase_compress_file .= "_dup";
    }
    if (lc $bad_snps_cnvs eq "y" or lc $bad_snps_cnvs eq "yes")
    {
        $ase_compress_file .= "_bsc";
    }
    if (lc $overlap eq "y" or lc $overlap eq "yes")
    {
        $ase_compress_file .= "_overlap";
    }
    
    $ase_compress_file .= ".tar.gz";
    
    chdir "$RNA_Path/$cds_plink";
    my @impute_hap_list = `ls All_Hets_Haps.chr*.sorted.bed TCGA*.snplist`;
    chdir $ase;
    
    open(IHO,">$impute_haps_archive");
    
    foreach my $impute_hap_files (@impute_hap_list)
    {
        print IHO"$RNA_Path\t$cds_plink/$impute_hap_files";
    }
    close (IHO);
    
    print "Making archive $impute_compress_file.\n";
    $parsing->archive_files("$impute_haps_archive","$impute_compress_file");
    `mv $RNA_Path/$impute_compress_file $Analysispath/$cancer_type/$finished_RNA`;
    chdir $ase;
    print "Making archive $ase_compress_file.\n";
    
    #sets up archive files before running archive_files routine
    open(MPA,">$RNA_CDS_Mpileups_archive");
    open(CDA,">$cds_sorted_archive");
    
    my @cds_mp;
    if (lc $overlap eq "y" or lc $overlap eq "yes")
    {
        @cds_mp =  `cat $RNA_CDS_Mpileups_overlap`;
        #get the paths to $ase_counts files
        $ase_analysis->print_ase_gene_level_path("$ase","$ase_counts","ase_counts_look.txt","ase_counts_getpath.txt","$ase_counts_archive",$overlap,$RNA_CDS_Mpileups_overlap);
        #get the paths to $ase_counts files
        $ase_analysis->print_ase_gene_level_path("$ase","$gene_level","gene_level_look.txt","gene_level_getpath.txt","$gene_level_archive",$overlap,$RNA_CDS_Mpileups_overlap);
    }
    else
    {
        @cds_mp =  `cat $RNA_CDS_Mpileups`;
        #get the paths to $ase_counts files
        $ase_analysis->print_ase_gene_level_path("$ase","$ase_counts","ase_counts_look.txt","ase_counts_getpath.txt","$ase_counts_archive",$overlap,$RNA_CDS_Mpileups);
        #get the paths to $ase_counts files
        $ase_analysis->print_ase_gene_level_path("$ase","$gene_level","gene_level_look.txt","gene_level_getpath.txt","$gene_level_archive",$overlap,$RNA_CDS_Mpileups);
    }
    
    foreach my $cd_mp (@cds_mp)
    {
        chomp ($cd_mp);
        my ($mpile,$tcga,$cds_bed) = split("\t",$cd_mp);
        print MPA "$ase\t$mpileups_dir/$mpile\n";
        print CDA "$RNA_Path\t$cds_sorted/$tcga.bed\n";
    }
    
    chdir $RNA_Path;
    my @Imputation_files = `ls $cancer_type\_TN_TCGA_Imputation.* $cancer_type\_TN_TCGA_All.*`;
    my @map_files = `ls $maps`;
    my @ped_files = `ls $peds`;
    my @phased_files = `ls $phased`;
    chdir $ase;
    
    open(IMO,">$Impute_archive");
    
    foreach my $impute_file(@Imputation_files)
    {
        print IMO "$RNA_Path\t$impute_file";
    }
    close (IMO);
    
    open(MO,">$maps_archive");
    
    foreach my $map_list (@map_files)
    {
        print MO "$RNA_Path\t$maps/$map_list";
    }
    close(MO);
    
    open(PO,">$peds_archive");
    
    foreach my $ped_list (@ped_files)
    {
        print PO "$RNA_Path\t$peds/$ped_list";
    }
    close(PO);
    
    open(PHO,">$phased_archive");
    
    foreach my $phased_list(@phased_files)
    {
        print PHO "$RNA_Path\t$phased/$phased_list";
    }
    close(PHO);
    
    my @matrix_list = `ls $ase/$matrix`;
    open(MO,">$matrix_archive");
    
    foreach my $matrix_files (@matrix_list)
    {
        print MO "$ase\t$matrix/$matrix_files";
    }
    close (MO);

    if (lc $bad_snps_cnvs eq "y" or lc $bad_snps_cnvs eq "yes")
    {
        $parsing->check_directory_existence("$ase/$problem_snps","$ase/$problem_cnvs");
        #sets up the $problem_snps and $problem_cnvs archives by using the $gene_level_archive file and substituting the $gene_level string with $problem_cnvs/$problem_snps and writing to the $problem_snps_archive and $problem_snps_archive files
        open(PCO,">$problem_cnv_archive");
        open(PSO,">$problem_snp_archive");
        my @gl_list = `cat $gene_level_archive`;
        
        foreach my $gl_files (@gl_list)
        {
            chomp($gl_files);
            my ($file_path,$sub_file) = split("\t",$gl_files);
            $sub_file =~ s/$gene_level/$problem_snps/;
            if (-e "$file_path/$sub_file")
            {
                print PSO "$file_path\t$sub_file\n";
                $sub_file =~ s/$problem_snps/$problem_cnvs/;
                print PCO "$file_path\t$sub_file\n";
            }
        }
        close(PSO);
        close(PCO);
        
        open(BSO,">$bad_snp_archive");
        open(BCO,">$bad_cnv_archive");
        my @cds = `cat $ase/$cds_sorted_archive`;
        
        foreach my $cd_files (@cds)
        {
            chomp($cd_files);
            my ($path,$filename) = split("\t",$cd_files);
            $filename =~ s/$cds_sorted/$bad_snps_bed/;
            if (-e "$path/$filename")
            {
                print BSO "$path\t$filename\n";
                $filename =~ s/$bad_snps_bed/$bad_cnvs/;
                print BCO "$path\t$filename\n";
            }
        }
        close (BSO);
        close (BCO);
        
        $parsing->archive_files("$ase/$maps_archive","$ase/$peds_archive","$ase/$phased_archive","$ase/$ase_counts_archive","$ase/$gene_level_archive","$ase/$RNA_CDS_Mpileups_archive","$ase/$cds_sorted_archive","$ase/$bad_snp_archive","$ase/$bad_cnv_archive","$ase/$problem_snp_archive","$ase/$problem_cnv_archive","$ase/$matrix_archive","$ase/$Impute_archive","$ase_compress_file");
    }
    else
    {
        $parsing->archive_files("$ase/$ase_counts_archive","$ase/$gene_level_archive","$ase/$RNA_CDS_Mpileups_archive","$ase/$cds_sorted_archive","$ase/$matrix_archive","$ase/$Impute_archive","$ase_compress_file");
    }
    
    mkdir "$Analysispath/$cancer_type/$finished_RNA" unless(-d "$Analysispath/$cancer_type/$finished_RNA");
    
    `mv $RNA_Path/$ase_compress_file $Analysispath/$cancer_type/$finished_RNA`;
}

print "All jobs have finished for $cancer_type.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
