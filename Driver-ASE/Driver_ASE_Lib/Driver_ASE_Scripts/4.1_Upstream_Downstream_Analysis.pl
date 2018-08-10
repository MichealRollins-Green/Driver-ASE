#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$Bin/..";
use WGS_Analysis;
use Parsing_Routines;
use Cwd "realpath";
use Getopt::Long;
use MCE::Map;
use strict;
use autodie;
no warnings 'once';

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $parsing = Driver_ASE_Lib::Parsing_Routines->new;
my $wgs_analysis = Driver_ASE_Lib::WGS_Analysis->new;

GetOptions(
    'cancer|c=s' => \my $cancer_type,#e.g. OV
    'overlap|o=s' => \my $overlap_data, #overlap option. only enter in if option was y|yes in script 4.0
    'remfiles|r=s' => \my $rem_files, #option to remove files
    'archive|a=s' => \my $archive, #option to archive data
    'overlapSelect|S=s' => \my $overlapSelect,#path to overlapSelect
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("4.1");

if ($help)
{
    $parsing->usage("4.1");
}

if (!defined $cancer_type)
{
    print STDERR "Cancer type was not entered!\n";
    $parsing->usage("4.1");
}

#default to no if option to overlap ids was not specified in the command line
if (!defined $overlap_data)
{
    $overlap_data = "no";
}
elsif (lc $overlap_data ne "y" and lc $overlap_data ne "yes" and lc $overlap_data ne "n" and lc $overlap_data ne "no")
{
    print STDERR "$overlap_data is an incorrect option! Must be yes|y|no|n.\n";
    $parsing->usage("4.1");
}

if (!defined $rem_files)
{
    $rem_files = "no";
}
elsif (lc $rem_files ne "y" and lc $rem_files ne "yes" and lc $rem_files ne "n" and lc $rem_files ne "no")
{
    print STDERR "$rem_files is an incorrect option! Must be yes|y|no|n.\n";
    $parsing->usage("4.1");   
}

if (!defined $archive)
{
    $archive = "no";
}
elsif (lc $archive ne "y" and lc $archive ne "yes" and lc $archive ne "n" and lc $archive ne "no")
{
    print STDERR "$archive is an incorrect option! Must be yes|y|no|n.\n";
    $parsing->usage("4.1");   
}

if (!defined $overlapSelect)
{
    $overlapSelect = "overlapSelect";
}

$parsing->CheckSoftware("$overlapSelect");

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $WGS_Path = "$Analysispath/$cancer_type/WGS_Analysis";
my ($wgs_dwnlds,$wgs_mpileups,$somatic,$annotate_vars) = ("wgs_dwnlds","wgs_mpileups","somatic_variants","annotate_vars");
my $finished_WGS = "$cancer_type\_finished_analysis_WGS";
my $overlap = "overlap";
my $reg = "reg";
my $tables = "$cancer_type\_tables";
my $somatic_calls = "somatic_calls";
my $mutations = "mutations";
my $annotations = "annotations";
my ($WGS_table_file,$WGS_table_overlap) = ("final_downloadtable_$cancer_type\_WGS.txt","final_downloadtable_$cancer_type\_WGS_overlap.txt");
my ($refseq_hg9_bed,$vars_bed,$ccds_bed,$ccds_tab,$stream_data) = ("refseq.ucsc.ensembl.mrna.hg9.nr.bed","vars.bed","ccds.hg19.bed","ccds.hg19.tab","stream_data");
my ($wgs_mpileup_archive,$somatic_archive,$annotations_archive,$somatic_calls_archive,$mutations_archive) = ("$wgs_mpileups\_varscan_to_archive.txt","$somatic\_to_archive.txt","$annotations\_to_archive.txt","$somatic_calls\_to_archive.txt","$mutations\_to_archive.txt");
my ($WGS_Mpileups_Full,$WGS_Mpileups_Overlapped) = ("WGS_Mpileups_Full.txt","WGS_Mpileups_Overlap.txt");
my $WGS_compress_file = "$cancer_type\_WGS_Analysis";

$parsing->check_directory_existence("$database_path","$Analysispath","$Analysispath/$cancer_type","$Analysispath/$cancer_type/$tables","$Analysispath/$cancer_type/$tables/$WGS_table_file","$WGS_Path","$WGS_Path/$somatic"); #check if directories or files exist
$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

chdir "$WGS_Path";

`mkdir -p $WGS_Path/$annotate_vars` unless(-d "$WGS_Path/$annotate_vars");
`find $WGS_Path/$annotate_vars/* 2>/dev/null|xargs rm -rf`;
#`rm -rf $WGS_Path/$annotate_vars/*`;

if (lc $overlap_data eq "y" or lc $overlap_data eq "yes")
{
    $parsing->check_directory_existence("$Analysispath/$cancer_type/$tables/$WGS_table_overlap");
    open (SOMO,">$WGS_Path/somatic_look.txt");
    my @get_som = `ls $WGS_Path/$somatic`;
    for my $som (@get_som)
    {
        chomp($som);
        my @sp_som = split("\\.",$som);
        print SOMO $sp_som[0],"\t",$sp_som[1],"\n";
    }
    close (SOMO);

    $parsing->vlookup("$WGS_Path/somatic_look.txt",1,"$Analysispath/$cancer_type/$tables/$WGS_table_overlap",1,1,"y","$WGS_Path/somatic_overlap_nan.txt");

    my @som_nan = `cat $WGS_Path/somatic_overlap_nan.txt`;
    open (SOML,">$WGS_Path/somatic_overlap_list.txt");
    foreach my $somo (@som_nan)
    {
        chomp($somo);
        if ($somo !~ /\tNaN+$/)
        {
            my @sp_som = split("\t",$somo);
            my $som_file_path = join(".",$sp_som[0],$sp_som[1]);
            print SOML "$WGS_Path/$somatic/$som_file_path","\n";
        }
    }
   close (SOML);

            #mk_files_for_wgs($WGS_Path/$somatic (directory that holds somatic variant files), output file)
   $wgs_analysis->mk_files_for_wgs("$WGS_Path/somatic_overlap_list.txt","$WGS_Path/somatic_list");
}
else
{
    #mk_files_for_wgs($WGS_Path/$somatic (directory that holds somatic variant files), output file)
    $wgs_analysis->mk_files_for_wgs("$WGS_Path/$somatic","$WGS_Path/somatic_list"); 
}

$parsing->matricize("$WGS_Path/somatic_list","$WGS_Path/somatic_list",1,6,"$WGS_Path");

mkdir "$WGS_Path/$mutations" unless (-d "$WGS_Path/$mutations");

`cp matrix.tab rowlabels.txt collabels.txt $WGS_Path/$mutations`;

#var_ID_to_bed(rowlabels.txt file created from matricize, output file)
$wgs_analysis->var_ID_to_bed("$WGS_Path/rowlabels.txt","$WGS_Path/$annotate_vars/$vars_bed");

chdir "$WGS_Path/$annotate_vars";

mkdir "$WGS_Path/$annotate_vars/$overlap" unless(-d "$WGS_Path/$annotate_vars/$overlap");
opendir (DD,"$database_path/$reg");
my @zcat = readdir (DD);
closedir (DD);
@zcat = grep{!/\.$/}@zcat;
my $bed;
mce_map
{
   $bed = "$_";
   $bed =~ s/\.gz//;
   #simpl(files from the $reg directory in the $database_path directory,output file)
   $wgs_analysis->simpl("zcat $database_path/$reg/$_ |","$WGS_Path/$annotate_vars/$overlap/$bed");
}@zcat;

#up_down_tss(the $refseq_hg9_bed file from the $database_path directory, upstream_of_tss_len, downstream_of_tss_len, output file)
$wgs_analysis->up_down_tss("$database_path/$refseq_hg9_bed",100,100,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_100");
#print1(output file created from up_down_tss, output file directed to the overlap directory)
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_100","$WGS_Path/$annotate_vars/$overlap/100.100.bed");

$wgs_analysis->up_down_tss("$database_path/$refseq_hg9_bed",500,500,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_500");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_500","$WGS_Path/$annotate_vars/$overlap/500.500.bed");

$wgs_analysis->up_down_tss("$database_path/$refseq_hg9_bed",1000,1000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_1kb");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_1kb","$WGS_Path/$annotate_vars/$overlap/1kb.1kb.bed");

$wgs_analysis->up_down_tss("$database_path/$refseq_hg9_bed",5000,1000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_5kb");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_5kb","$WGS_Path/$annotate_vars/$overlap/5kb.1kb.bed");

$wgs_analysis->up_down_tss("$database_path/$refseq_hg9_bed",10000,10000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10kb");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10kb","$WGS_Path/$annotate_vars/$overlap/10kb.1kb.bed");

$wgs_analysis->up_down_tss("$database_path/$refseq_hg9_bed",50000,1000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_50kb");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_50kb","$WGS_Path/$annotate_vars/$overlap/50kb.1kb.bed");

$wgs_analysis->up_tss_gene("$database_path/$refseq_hg9_bed",10000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000","$WGS_Path/$annotate_vars/$overlap/10kb.gene.bed");

$parsing->pull_column("$WGS_Path/$annotate_vars/$vars_bed",4,"$WGS_Path/$annotate_vars/$stream_data"."1");

#vlookem_all($WGS_Path/$annotate_vars/$overlap (directory that stores files created from simpl routine, $WGS_Path/$annotate_vars/$vars_bed (path to $vars_bed file), $WGS_Path/$annotate_vars (path to directory that holds annotation data), $overlapSelect (overlapSelect command), $stream_data (name of input file))
$wgs_analysis->vlookem_all("$WGS_Path/$annotate_vars/$overlap","$WGS_Path/$annotate_vars/$vars_bed","$WGS_Path/$annotate_vars","$overlapSelect","$stream_data");

opendir (DD,"$WGS_Path/$annotate_vars");
my @dd = readdir (DD);
@dd = grep{/^$stream_data*/}@dd;
closedir (DD);
my $sdata_num;
my $highest_sdata_num = 0;
my @a;
for (my $i = 0;$i < scalar(@dd);$i++)
{
   @a = split (/(\d+)/,$dd[$i]);
   $sdata_num = $a[1];
   if ($highest_sdata_num < $sdata_num)
   { 
       $highest_sdata_num = $sdata_num;
   }
}

my $sdata = $a[0] . $highest_sdata_num;

$parsing->pull_column("$WGS_Path/$annotate_vars/$sdata",1,"$WGS_Path/$annotate_vars/rowlabels.txt");
$parsing->pull_column("$WGS_Path/$annotate_vars/$sdata","2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19","$WGS_Path/$annotate_vars/matrix.tab");

#my_up_down_tss($database_path/$refseq_hg9_bed ($refseq_hg9_bed file from the $database_path directory), upstream_of_tss_len, downstream_of_tss_len, output file)
$wgs_analysis->my_up_down_tss("$database_path/$refseq_hg9_bed",10000,0,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000_0_pull");
$parsing->pull_column("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000_0_pull","1,2,3,4","$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000_0_sort");
`sort -k 1,1 -k 2,2n $WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000_0_sort > $WGS_Path/$annotate_vars/gene.bed`;
`$overlapSelect $WGS_Path/$annotate_vars/$vars_bed $WGS_Path/$annotate_vars/gene.bed -idOutput $WGS_Path/$annotate_vars/vars_gene_out`;

$parsing->strip_head("$WGS_Path/$annotate_vars/vars_gene_out","$WGS_Path/$annotate_vars/vars_gene_out_pull");
$parsing->pull_column("$WGS_Path/$annotate_vars/vars_gene_out_pull","2,1","$WGS_Path/$annotate_vars/vars_gene_out_sort");
`sort -k 1,1 $WGS_Path/$annotate_vars/vars_gene_out_sort > $WGS_Path/$annotate_vars/var2gene.tab`;

mkdir "$WGS_Path/$annotate_vars/$somatic_calls" unless (-d "$WGS_Path/$annotate_vars/$somatic_calls");
`cp $WGS_Path/$annotate_vars/var2gene.tab $WGS_Path/$annotate_vars/$somatic_calls`;

#mk_coding($WGS_Path/$annotate_vars/$vars_bed ($vars_bed file), output file)
$wgs_analysis->mk_coding("$WGS_Path/$annotate_vars/$vars_bed","$WGS_Path/$annotate_vars/vars_grep");
`cat $WGS_Path/$annotate_vars/vars_grep | grep shift > $WGS_Path/$annotate_vars/vars_sort.txt`;
`sort -k 1,1 -k 2,2n $WGS_Path/$annotate_vars/vars_sort.txt > $WGS_Path/$annotate_vars/shifts.bed`;
$wgs_analysis->mk_coding("$WGS_Path/$annotate_vars/$vars_bed","$WGS_Path/$annotate_vars/vars_grep");
`cat $WGS_Path/$annotate_vars/vars_grep | grep -v shift > $WGS_Path/$annotate_vars/vars_sort.txt`;
`sort -k 1,1 -k 2,2n $WGS_Path/$annotate_vars/vars_sort.txt > $WGS_Path/$annotate_vars/subs.bed`;

#Syn($WGS_Path/$annotate_vars/subs.bed (subs bed file), genome, output file , path to the $database_path directory, $overlapSelect (overlapSelect command), $ccds_bed (ccds bed file in $database_path directory), $ccds_tab (ccds tab file in $database_path directory))
$wgs_analysis->Syn("$WGS_Path/$annotate_vars/subs.bed","hg19","$WGS_Path/$annotate_vars/subs.syn","$database_path","$overlapSelect","$ccds_bed","$ccds_tab");
#mk_uniq_label(output file from Syn,output file)
$wgs_analysis->mk_uniq_label("$WGS_Path/$annotate_vars/subs.syn","$WGS_Path/$annotate_vars/subs.syn_sort");
`sort -k 1,1 -k 12,12n $WGS_Path/$annotate_vars/subs.syn_sort > $WGS_Path/$annotate_vars/subs.syn_sort2`;
`sort -k 1,1 -u $WGS_Path/$annotate_vars/subs.syn_sort2 > $WGS_Path/$annotate_vars/subs.syn2`;

$parsing->pull_column("$WGS_Path/$annotate_vars/subs.syn2","2,3,4,5,6,7,8,9,10,11,12,13,14","$WGS_Path/$annotate_vars/subs.syn2_process");

#process_syn(output file from pull_column,output file)
$wgs_analysis->process_syn("$WGS_Path/$annotate_vars/subs.syn2_process","$WGS_Path/$annotate_vars/subs.syn2_processed");

`$overlapSelect $database_path/$ccds_bed $WGS_Path/$annotate_vars/shifts.bed $WGS_Path/$annotate_vars/overlap_output`;

#process_shifts(output file from overlapSelect,output file)
$wgs_analysis->process_shifts("$WGS_Path/$annotate_vars/overlap_output","$WGS_Path/$annotate_vars/subs.syn2_processed");

#pt(the rowlabels.txt file created from matricize in script 4.1,output file)
$wgs_analysis->pt("$WGS_Path/rowlabels.txt","$WGS_Path/$annotate_vars/rowlabels_processed.txt");

$parsing->vlookup("$WGS_Path/$annotate_vars/rowlabels_processed.txt",1,"$WGS_Path/$annotate_vars/subs.syn2_processed",1,"2,3","y","$WGS_Path/$annotate_vars/syn_pull");

$parsing->pull_column("$WGS_Path/$annotate_vars/syn_pull","2,3,4","$WGS_Path/$annotate_vars/syn.tab");

mkdir "$WGS_Path/$annotate_vars/$annotations" unless (-d "$WGS_Path/$annotate_vars/$annotations");
`cp matrix.tab rowlabels.txt collabels.txt syn.tab $WGS_Path/$annotate_vars/$annotations`;

`$overlapSelect $database_path/$ccds_bed $WGS_Path/$annotate_vars/$vars_bed $WGS_Path/$annotate_vars/vars_ccds.bed`;
`$overlapSelect $WGS_Path/$annotate_vars/vars_ccds.bed $database_path/$refseq_hg9_bed -idOutput $WGS_Path/$annotate_vars/vars_refseq_bed`;

$parsing->strip_head("$WGS_Path/$annotate_vars/vars_refseq_bed","$WGS_Path/$annotate_vars/vars_refseq_bed_pull");
$parsing->pull_column("$WGS_Path/$annotate_vars/vars_refseq_bed_pull","2,1","$WGS_Path/$annotate_vars/var2gene.coding.tab");

`cp $WGS_Path/$annotate_vars/var2gene.coding.tab $WGS_Path/$annotate_vars/$somatic_calls`;

if (lc $overlap_data eq "y" or lc $overlap_data eq "yes")
{
    $parsing->pull_column("$Analysispath/$cancer_type/$tables/$WGS_table_overlap","1,2","$WGS_Path/$annotate_vars/look.txt");
}
else
{
    $parsing->pull_column("$Analysispath/$cancer_type/$tables/$WGS_table_file","1,2","$WGS_Path/$annotate_vars/look.txt");     
}

#join UUID and TCGA ID and print the fields
my @UUID_TCGA =  `cat look.txt`;
open(LT,">look.tab");
foreach my $ut (@UUID_TCGA)
{
    chomp($ut);
    my @sp_ut = split("\t",$ut);
    my $utn = join(".",@sp_ut);
    print LT "$ut\t$utn\n";
}
close (LT);

`cp $WGS_Path/$annotate_vars/look.tab $WGS_Path/$annotate_vars/$somatic_calls`;

#two files: mut_ase_look.txt and coding_genes.tab are pre-created
`cp $database_path/mut_ase_look.txt $WGS_Path/$annotate_vars/$somatic_calls`;
`cp $database_path/coding_genes.tab $WGS_Path/$annotate_vars/$somatic_calls`;

if (lc $rem_files eq "y" or lc $rem_files eq "yes")
{
    chdir "$WGS_Path";
    my @del_files = $parsing->get_only_files_in_dir("$WGS_Path");
    #Removes all files after the processes have finished.
    for (my $i = 0;$i < scalar(@del_files);$i++)
    {
        `rm "$del_files[$i]"`;
    }
    undef @del_files;
    chdir "$WGS_Path/$annotate_vars";
    @del_files = $parsing->get_only_files_in_dir("$WGS_Path/$annotate_vars");
    for (my $i = 0;$i < scalar(@del_files);$i++)
    {
        `rm "$del_files[$i]"`;
    }
}

if (lc $archive eq "y" or lc $archive eq "yes")
{
    chdir $WGS_Path;
    my @wgs_mps;
    if (lc $overlap_data eq "y" or lc $overlap_data eq "yes")
    {
        $WGS_compress_file .= "_overlap.tar.gz";
       open (WMO,">$WGS_Path/$wgs_mpileup_archive");
       @wgs_mps = `cat $WGS_Path/$WGS_Mpileups_Overlapped`;
    }
    else
    {
        open (WMO,">$WGS_Path/$wgs_mpileup_archive");
        $WGS_compress_file .= ".tar.gz";
        @wgs_mps = `cat $WGS_Path/$WGS_Mpileups_Full`;
    }
    
    print "Making archive $WGS_compress_file.\n";
    
    open (WSO,">$somatic_archive");
    my $prev_som = "";
    foreach my $wgs_lists (@wgs_mps)
    {
        chomp($wgs_lists);
        my ($mpileups,$somatic_files,$file_path,$wgs_filename) = split("\t",$wgs_lists);
        print WMO "$file_path\t$wgs_filename\n";
        $file_path =~ s/$file_path/$WGS_Path/;
        #Only print if the currnet file name is not the same as the previous one
        if ($somatic_files ne $prev_som)
        {
            print WSO "$file_path\t$somatic/$somatic_files\n";
            $prev_som = $somatic_files;
        }
    }
    close (WMO);
    close (WSO);
    
    chdir "$WGS_Path";
    
    my @mutation = `ls $mutations`;
    open (MUTO,">$WGS_Path/$mutations_archive");
    
    foreach my $mut (@mutation)
    {
        print MUTO "$WGS_Path\t$mutations/$mut";
    }
    close (MUTO);
    
    chdir "$WGS_Path/$annotate_vars";
    
    open (ANNO,">$WGS_Path/$annotations_archive");
    my @annotate_files = `ls $annotations`;
    
    foreach my $ann_list (@annotate_files)
    {
        print ANNO "$WGS_Path/$annotate_vars\t$annotations/$ann_list";
    }
    close (ANNO);
    
    open (SOMC,">$WGS_Path/$somatic_calls_archive");
    
    my @som_calls = `ls $somatic_calls`;
    
    foreach my $sc (@som_calls)
    {
        print SOMC "$WGS_Path/$annotate_vars\t$somatic_calls/$sc";
    }
    close (SOMC);
    
    chdir $WGS_Path;
    $parsing->archive_files("$WGS_Path/$wgs_mpileup_archive","$WGS_Path/$somatic_archive","$WGS_Path/$annotations_archive","$WGS_Path/$mutations_archive","$WGS_Path/$somatic_calls_archive","$WGS_compress_file");
    
    mkdir "$Analysispath/$cancer_type/$finished_WGS" unless (-d "$Analysispath/$cancer_type/$finished_WGS");

    `mv $WGS_Path/$annotate_vars/$WGS_compress_file $Analysispath/$cancer_type/$finished_WGS`;
}

print "All jobs have finished for $cancer_type.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
