#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$Bin/..";
use Parsing_Routines;
use WGS_Analysis;
use Cwd "realpath";
use Getopt::Long;
use MCE::Map;
use strict;

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $parsing = TCGA_Lib::Parsing_Routines->new;
my $wgs_analysis = TCGA_Lib::WGS_Analysis->new;

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage;

if($help)
{
    $parsing->usage;
}

my $TCGA_Pipeline_Dir = realpath("../../");
my $database_path = "$TCGA_Pipeline_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $WGS_Path = "$Analysispath/$disease_abbr/WGS_Analysis";
my $somatic = "somatic_variants";
my $annotate_vars = "annotate_vars";
my $finished_WGS = "$disease_abbr\_finished_analysis_WGS";
my $overlap = "overlap";
my $reg = "reg";

if(!defined $disease_abbr)
{
    print "disease type was not entered!\n";
    $parsing->usage;
}

#Check if there is not Database directory
if(!(-d "$database_path"))
{
    print STDERR "$database_path does not exist, it was either moved, renamed, deleted or has not been downloaded.\nPlease check the README.md file on the github page to find out where to get the Database directory.\n";
    exit;
}

#Checks if there is not Analysis directory.
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

if(-d "$Analysispath/$disease_abbr/$disease_abbr\_tables")
{
   if (!(-e "$Analysispath/$disease_abbr/$disease_abbr\_tables/final_downloadtable_$disease_abbr\_WGS.txt"))
    {
        print STDERR "$Analysispath/$disease_abbr/$disease_abbr\_tables/final_downloadtable_$disease_abbr\_WGS.txt does not exist. It was either moved, renamed or deleted.\n";
        print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
        exit;
    } 
}
else
{
    print STDERR "$Analysispath/$disease_abbr/$disease_abbr\_tables does not exist. It was either moved, renamed or deleted,\n";
    print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
    exit;
}

if(!(-d $WGS_Path))
{
    print STDERR "$WGS_Path does not exits, it was either deleted, moved or renamed.\n";
    print STDERR "Please run script 4.0_Somatic_Variants.pl.\n";
    exit;
}

chdir "$WGS_Path";

if (!(-d "$WGS_Path/$somatic"))
{
    print STDERR "$WGS_Path/$somatic does not exist. It was either deleted, moved or renamed.\n";
    print STDERR "Please run script 4.0_Somatic_Variants.pl.\n";
    exit;
}

`mkdir -p $WGS_Path/$annotate_vars` unless(-d "$WGS_Path/$annotate_vars");

#mk_files_for_wgs(somatic_variants directory,output file)
$wgs_analysis->mk_files_for_wgs("$WGS_Path/$somatic","$WGS_Path/ff");

$parsing->matricize("$WGS_Path/ff","$WGS_Path/ff",1,6,"$WGS_Path");

`tar -zcvf mutations.tar.gz matrix.tab rowlabels.txt collabels.txt`;

mkdir "$Analysispath/$disease_abbr/$finished_WGS" unless(-d "$finished_WGS");
`mv $WGS_Path/mutations.tar.gz $Analysispath/$disease_abbr/$finished_WGS`;

#var_ID_to_bed(rowlabels.txt file created from matricize,output file)
$wgs_analysis->var_ID_to_bed("$WGS_Path/rowlabels.txt","$WGS_Path/$annotate_vars/vars.bed");

chdir "$WGS_Path/$annotate_vars";

mkdir "$WGS_Path/$annotate_vars/$overlap";
opendir(DD,"$database_path/$reg") or die "Can't open $database_path/$reg: $!\n";
my @zcat = readdir(DD);
closedir(DD);
@zcat = grep{!/\.$/}@zcat;
my $bed;
mce_map
{
    $bed = "$_";
    $bed =~ s/\.gz//;
    #simpl(files from the reg directory in the Database directory,output file)
    $wgs_analysis->simpl("zcat $database_path/$reg/$_ |","$WGS_Path/$annotate_vars/$overlap/$bed");
}@zcat;

#up_down_tss(the refseq.ucsc.ensembl.mrna.hg9.nr.bed file from the Database directory,upstream_of_tss_len,downstream_of_tss_len,output file)
$wgs_analysis->up_down_tss("$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",100,100,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_100");
#print1(output file created from up_down_tss,output file directed to the overlap directory)
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_100","$WGS_Path/$annotate_vars/$overlap/100.100.bed");

$wgs_analysis->up_down_tss("$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",500,500,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_500");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_500","$WGS_Path/$annotate_vars/$overlap/500.500.bed");

$wgs_analysis->up_down_tss("$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",1000,1000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_1kb");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_1kb","$WGS_Path/$annotate_vars/$overlap/1kb.1kb.bed");

$wgs_analysis->up_down_tss("$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",5000,1000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_5kb");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_5kb","$WGS_Path/$annotate_vars/$overlap/5kb.1kb.bed");

$wgs_analysis->up_down_tss("$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",10000,10000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10kb");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10kb","$WGS_Path/$annotate_vars/$overlap/10kb.1kb.bed");

$wgs_analysis->up_down_tss("$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",50000,1000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_50kb");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_50kb","$WGS_Path/$annotate_vars/$overlap/50kb.1kb.bed");

$wgs_analysis->up_tss_gene("$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",10000,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000");
$wgs_analysis->print_1("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000","$WGS_Path/$annotate_vars/$overlap/10kb.gene.bed");

$parsing->pull_column("$WGS_Path/$annotate_vars/vars.bed",4,"$WGS_Path/$annotate_vars/tt1");

#vlookem_all(input files from the overlap directory,vars.bed file,path to the user defined directory from command line or default)
$wgs_analysis->vlookem_all("$WGS_Path/$annotate_vars/$overlap","$WGS_Path/$annotate_vars/vars.bed","$WGS_Path/$annotate_vars");

opendir(DD,"$WGS_Path/$annotate_vars") or die "Can't open $WGS_Path/$annotate_vars: $!\n";
my @dd = readdir(DD);
@dd = grep{/^tt*/}@dd;
closedir(DD);
my $tt_num;
my $highest_tt = 0;
my @a;
for(my $i = 0;$i < scalar(@dd);$i++)
{
    @a = unpack("(a2)*",$dd[$i]);
    $tt_num = $a[1];
    if($highest_tt < $tt_num)
    { 
        $highest_tt = $tt_num;
    }
}

my $tt = $a[0] . $highest_tt;

$parsing->pull_column("$WGS_Path/$annotate_vars/$tt",1,"$WGS_Path/$annotate_vars/rowlabels.txt");
$parsing->pull_column("$WGS_Path/$annotate_vars/$tt","2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19","$WGS_Path/$annotate_vars/matrix.tab");

#my_up_down_tss(refseq.ucsc.ensembl.mrna.hg9.nr.bed file from the Database directory,upstream_of_tss_len,downstream_of_tss_len,output file)
$wgs_analysis->my_up_down_tss("$database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed",10000,0,"$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000_0_pull");
$parsing->pull_column("$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000_0_pull","1,2,3,4","$WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000_0_sort");
`sort -k 1,1 -k 2,2n $WGS_Path/$annotate_vars/mrna.hg9.nr.bed_10000_0_sort > $WGS_Path/$annotate_vars/gene.bed`;
`overlapSelect $WGS_Path/$annotate_vars/vars.bed $WGS_Path/$annotate_vars/gene.bed -idOutput $WGS_Path/$annotate_vars/out`;

$parsing->strip_head("$WGS_Path/$annotate_vars/out","$WGS_Path/$annotate_vars/out_pull");
$parsing->pull_column("$WGS_Path/$annotate_vars/out_pull","2,1","$WGS_Path/$annotate_vars/out_sort");
`sort -k 1,1 $WGS_Path/$annotate_vars/out_sort > $WGS_Path/$annotate_vars/var2gene.tab`;

#mk_coding(vars.bed file,output file)
$wgs_analysis->mk_coding("$WGS_Path/$annotate_vars/vars.bed","$WGS_Path/$annotate_vars/vars_grep");
`cat $WGS_Path/$annotate_vars/vars_grep | grep shift > $WGS_Path/$annotate_vars/vars_sort.txt`;
`sort -k 1,1 -k 2,2n $WGS_Path/$annotate_vars/vars_sort.txt > $WGS_Path/$annotate_vars/shifts.bed`;
$wgs_analysis->mk_coding("$WGS_Path/$annotate_vars/vars.bed","$WGS_Path/$annotate_vars/vars_grep");
`cat $WGS_Path/$annotate_vars/vars_grep | grep -v shift > $WGS_Path/$annotate_vars/vars_sort.txt`;
`sort -k 1,1 -k 2,2n $WGS_Path/$annotate_vars/vars_sort.txt > $WGS_Path/$annotate_vars/subs.bed`;

#Syn(subs.ved file,genome,output file,path to the Database directory)
$wgs_analysis->Syn("$WGS_Path/$annotate_vars/subs.bed","hg19","$WGS_Path/$annotate_vars/subs.syn","$database_path");
#mk_uniq_label(output file from Syn,output file)
$wgs_analysis->mk_uniq_label("$WGS_Path/$annotate_vars/subs.syn","$WGS_Path/$annotate_vars/subs.syn_sort");
`sort -k 1,1 -k 12,12n $WGS_Path/$annotate_vars/subs.syn_sort > $WGS_Path/$annotate_vars/subs.syn_sort2`;
`sort -k 1,1 -u $WGS_Path/$annotate_vars/subs.syn_sort2 > $WGS_Path/$annotate_vars/subs.syn2`;

$parsing->pull_column("$WGS_Path/$annotate_vars/subs.syn2","2,3,4,5,6,7,8,9,10,11,12,13,14","$WGS_Path/$annotate_vars/subs.syn2_process");

#process_syn(output file from pull_column,output file)
$wgs_analysis->process_syn("$WGS_Path/$annotate_vars/subs.syn2_process","$WGS_Path/$annotate_vars/t1");

`overlapSelect $database_path/ccds.hg19.bed $WGS_Path/$annotate_vars/shifts.bed $WGS_Path/$annotate_vars/o`;

#process_shifts(output file from overlapSelect,output file)
$wgs_analysis->process_shifts("$WGS_Path/$annotate_vars/o","$WGS_Path/$annotate_vars/t1");

#pt(the rowlabels.txt file created from matricize in script 4.1,output file)
$wgs_analysis->pt("$WGS_Path/rowlabels.txt","$WGS_Path/$annotate_vars/t2");

$parsing->vlookup("$WGS_Path/$annotate_vars/t2",1,"$WGS_Path/$annotate_vars/t1",1,"2,3","y","$WGS_Path/$annotate_vars/syn_pull");

$parsing->pull_column("$WGS_Path/$annotate_vars/syn_pull","2,3,4","$WGS_Path/$annotate_vars/syn.tab");

`tar zcvf annotations.tar.gz matrix.tab rowlabels.txt collabels.txt syn.tab`;
mkdir "$Analysispath/$disease_abbr/$finished_WGS" unless(-d "$Analysispath/$disease_abbr/$finished_WGS");
`mv $WGS_Path/$annotate_vars/annotations.tar.gz $Analysispath/$disease_abbr/$finished_WGS`;

`overlapSelect $database_path/ccds.hg19.bed $WGS_Path/$annotate_vars/vars.bed $WGS_Path/$annotate_vars/o.bed`;
`overlapSelect $WGS_Path/$annotate_vars/o.bed $database_path/refseq.ucsc.ensembl.mrna.hg9.nr.bed -idOutput $WGS_Path/$annotate_vars/out`;

$parsing->strip_head("$WGS_Path/$annotate_vars/out","$WGS_Path/$annotate_vars/out_pull");
$parsing->pull_column("$WGS_Path/$annotate_vars/out_pull","2,1","$WGS_Path/$annotate_vars/var2gene.coding.tab");

$parsing->pull_column("$Analysispath/$disease_abbr/$disease_abbr\_tables/final_downloadtable_$disease_abbr\_WGS.txt","1,2","$WGS_Path/$annotate_vars/look.tab");

#two files: mut_ase_look.txt and coding_genes.tab are pre-created
`cp $database_path/mut_ase_look.txt $WGS_Path/$annotate_vars`;
`cp $database_path/coding_genes.tab $WGS_Path/$annotate_vars`;
`tar -zcvf somatic_calls.tar.gz look.tab var2gene.tab var2gene.coding.tab coding_genes.tab mut_ase_look.txt`;
`mv $WGS_Path/$annotate_vars/somatic_calls.tar.gz $Analysispath/$disease_abbr/$finished_WGS`;

chdir "$WGS_Path";
my @del_files = $parsing->get_only_files_in_dir("$WGS_Path");
#Removes all files after the processes have finished.
for(my $i = 0;$i < scalar(@del_files);$i++)
{
    `rm "$del_files[$i]"`;
}
undef @del_files;
chdir "$WGS_Path/$annotate_vars";
@del_files = $parsing->get_only_files_in_dir("$WGS_Path/$annotate_vars");
for(my $i = 0;$i < scalar(@del_files);$i++)
{
    `rm "$del_files[$i]"`;
}

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
