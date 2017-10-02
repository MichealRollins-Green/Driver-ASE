#!/usr/bin/perl
package Driver_ASE_Lib::Bad_SNPS_and_CNVS;

use FindBin qw($Bin);
use strict;
use warnings;
use lib "$Bin/";
use MCE::Map;
use Parsing_Routines;
use autodie;
no warnings 'once';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

my $parsing = Driver_ASE_Lib::Parsing_Routines->new;

###############################################################################
# Bad_SNPS_and_CNVS.pm version 1
# Written by Tomas Babak, Zhongshan Cheng, Michael Vermeulen and Micheal Rollins-Green
#  in the Department of Biology at Queen's University.
#
# Bad_SNPS_and_CNVS
#
# Module interface
#
###############################################################################

###############################################################################
# Constructor
###############################################################################
sub new 
{
    my $this = shift;
    
    my $class = ref($this) || $this;

    my $self = {};
    bless $self,$class;
    return($self);
}

##############################################subs#############################
sub pull_bed
{
    my $annot_csv_file = shift; #file that was created from GenomeWideSNP_6.cn.na35.annot.csv.zip
    my $self = $annot_csv_file and $annot_csv_file = shift if ref $annot_csv_file;
    my $outfile = shift; #file where data will be outputted
    
    open (PB,$annot_csv_file);
    open (PBO,">$outfile");
    while(my $r = <PB>)
    {
        if($r=~/^\#/)
        {
            next;
        }
        chomp($r);
        $r =~ s/\"//g;
        my @a = split(",",$r);
        print PBO "chr".$a[1], "\t", $a[2], "\t", $a[3], "\t", $a[0], "\n";
    }
    close (PB);
    close (PBO);
}

sub pull_normal
{
    my $cnv_file_list = shift; #file that has the list of the downloaded Copy number estimate files downloaded in 0_Download_SNPArray_From_GDC.pl script
    my $self = $cnv_file_list and $cnv_file_list = shift if ref $cnv_file_list;
    my $normal_cnv = shift; #outfile for normal cnvs
    
    open (PN,$cnv_file_list);
    open (PNO,">$normal_cnv");

    while(my $r = <PN>)
    {
        chomp($r);
        my $tcga = [split("\\.",$r)]->[-1];#Get TCGA ID
        my $tumor_type = [split("-",$tcga)]->[-1];#Get tumor/normal sample type
        $tumor_type =~ s/[a-zA-Z]$//;
        $tcga =~ s/-[0-9]+[a-zA-Z]$//;
        if($tumor_type > 9)
        {
            print PNO $r, "\t", $tcga, "\t", $tumor_type, "\n";
        }
    }
    close (PN);
    close (PNO);
}

sub mv_bcode
{
    my $cds_sorted_list = shift; #file that contains list of sorted cds bed files
    my $self = $cds_sorted_list and $cds_sorted_list = shift if ref $cds_sorted_list;
    my $cds_seperated = shift; #out file that breaks up the bed file name
    
    open (MVB, $cds_sorted_list);
    open (MVBO, ">$cds_seperated");
    while(my $r = <MVB>)
    {
        chomp($r);
        my $TCGA_ID = $r;
        my $sample = [split("-",$TCGA_ID)]->[-1];#get sample type
        $TCGA_ID =~ s/-[0-9]+[a-zA-Z].[a-zA-Z]+$//;#Get rid of sample and .bed extension
        $sample =~ s/[a-zA-Z].[a-zA-Z]+$//;#Get rid of sample letter and .bed extension
        print MVBO $r, "\t", $TCGA_ID, "\t", $sample, "\n";
    }
    close (MVB);
    close (MVBO);
}

sub hv_cd
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open (HVC,$infile);
    open (HVCO,">$outfile");
    while(my $r = <HVC>)
    {
        chomp($r);
        my @a = split("\t",$r);
        unless($a[3] eq 'NaN')
        {
            print HVCO $r, "\n";
        }
    }
    close (HVC);
    close (HVCO);
}

sub dump_non_01_10_11
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open (DNI,$infile);
    open (DNO,">$outfile");
    while(my $r = <DNI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        if($a[3] eq '06')
        {
            next;
        }
        print DNO $r, "\n";
    }
    close (DNI);
    close (DNO);
}

sub mk_tn_tables
{
    my $snp_raw_files_dir = shift; #directory where Genotype/Copy number estimate data was downloaded
    my $self = $snp_raw_files_dir and $snp_raw_files_dir = shift if ref $snp_raw_files_dir;
    my $raw_file_name = shift; #raw file name(same as directory name but will be used for files)
    
    `ls "$snp_raw_files_dir" > "$raw_file_name"\_raw_files`;

    open (HC,"$raw_file_name\_raw_files");
    open (HCT,">$raw_file_name\_T");
    open (HCN,">$raw_file_name\_N");
    
    while(my $r = <HC>)
    {
        chomp($r);
        my ($UUID,$TCGA) = split("\\.",$r);
        my $tcga_id = $TCGA;
        $tcga_id =~ s/-[0-9]+[a-zA-Z]$//;
        my $tumor_type = [split("-",$TCGA)]->[-1];
        $tumor_type =~ s/[a-zA-Z]$//;
        if ($tumor_type > 9)
        {
            print HCN $TCGA, "\t", $r, "\t", $tcga_id, "\t", $tumor_type, "\n"; 
        }
        elsif($tumor_type <= 9)
        {
            print HCT $TCGA, "\t", $r, "\t", $tcga_id, "\t", $tumor_type, "\n";   
        }
        else
        {
            print STDERR "Something is wrong with $tumor_type!\n";
            exit;
        }
    }
    close (HCN);
    close (HCT);
}

sub run_snps
{
    my $raw_snp_dir = shift; #Genotypes and Copy number estimate directory
    my $self = $raw_snp_dir and $raw_snp_dir = shift if ref $raw_snp_dir;
    my ($cnv_geno_tn_file,$RNA_Path) = @_; #file with tn cnvs/genos, path to RNA analysis data
    
    my ($bird,$copy) = split(",",$raw_snp_dir);
    $bird =~ s/\/$//;
    $copy =~ s/\/$//;
    mce_map_f
    {
        chomp($_);
        my @a = split("\t",$_);
        print "The job for $a[0] was submitted.\n";# ? $a[3] is TCGA ID;
        my $tmp = $parsing->temp_filename();
        if (-f "$bird/$a[1]" and -f "$bird/$a[2]")
        {
            $parsing->vlookup("$bird/$a[2]",1,"$bird/$a[1]",1,2,"y","$tmp.txt");
            flag_snps("$tmp.txt","$RNA_Path/$a[0]");
        }
        elsif(-f "$copy/$a[1]" and -f "$copy/$a[2]")
        {
            $parsing->vlookup("$copy/$a[2]",1,"$copy/$a[1]",1,2,"y","$tmp.txt");
            flag_snps("$tmp.txt","$RNA_Path/$a[0]"); 
        }
        else
        {
            print "SNP array data not available for $a[1].\n";
        }
        `rm -f $tmp.txt`;
    }$cnv_geno_tn_file; 
}

sub flag_snps
{
    my $infile = shift; #tmp file
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift; #path to output file
    
    open (FSI,$infile);
    open (FSO,">$outfile");
    #Remove two unuseful lines;
    my $r = <FSI>;
    $r = <FSI>;
    
    while($r = <FSI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        if($a[1] ne $a[3])
        {
            print FSO $a[0], "\n";
        } 
    }
    close (FSI);
    close (FSO);
}

sub get_cd_bed
{
    my $snp6_cd_file = shift; #Path to snp6.cd.txt
    my $self = $snp6_cd_file and $snp6_cd_file = shift if ref $snp6_cd_file;
    my ($bad_snps_bed,$bad_snps) = @_; #bad snps bed directory, bad snps directory
    
    opendir (DD,"$bad_snps");
    my @dd = readdir(DD);
    closedir (DD);
    @dd = grep{!/^\./ && -f "$bad_snps/$_"}@dd;#gets rid off . and ..
    mce_map
    {
        print "Working on $_\n";
        
        my $tmp = $parsing->temp_filename();
       `sort -u -k 1,1 $bad_snps/$_ > $bad_snps/$_.uniq`;
       $parsing->vlookup("$bad_snps/$_.uniq",1,"$snp6_cd_file",4,"1,2,3","y","$bad_snps/$tmp.1.txt");
       `grep NaN -v $bad_snps/$tmp.1.txt > $bad_snps/$tmp.2.txt`;
       $parsing->pull_column("$bad_snps/$tmp.2.txt","2,3,4","$bad_snps/$tmp.3.txt");
       `sort -k 1,1 -k 2,2n $bad_snps/$tmp.3.txt > $bad_snps_bed/$_`;
       `rm $bad_snps/$tmp.*.txt`;
    }@dd;
}

sub run_cnv
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my ($copy,$bird,$cnv_bed,$RNA_Path_cnvs,$temp_dir) = @_; #cnv directory, genotypes directory, cnv.hg19.bed file, path to directory where bad cnvs with be written to, temp directory

    open (my $CIN,"$infile");
    my @CNVS = <$CIN>;

    mce_map
    {
        chomp($_);
        my @a = split("\t",$_);
        print "The job for $a[0] was submitted.\n";
        
        my $rr = rand();
        $rr = substr($rr,2,6);
        #This will check if the CNV files exist. They should always exist unless there were changes made to the directory prior to this running.
        if (-e "$copy/$a[5]" and -e "$copy/$a[6]")
        {
            $parsing->vlookup("$copy/$a[5]",1,"$copy/$a[6]",1,2,"y","$temp_dir/t.$rr");
            $parsing->vlookup("$temp_dir/t.$rr",1,$cnv_bed,4,"1,2,3","y","$temp_dir/t.$rr.1.txt");
            $parsing->strip_head("$temp_dir/t.$rr.1.txt","$temp_dir/t.$rr.2.txt",2);
            `sort -k 4,4 -k 5,5n $temp_dir/t.$rr.2.txt > $temp_dir/t.$rr.3.txt`;
            smooth("$temp_dir/t.$rr.3.txt","$temp_dir/t.$rr.4.txt");
            delin_cnv("$temp_dir/t.$rr.4.txt","$temp_dir/t.$rr.5.txt",0.5);
            `grep NaN -v $temp_dir/t.$rr.5.txt > $RNA_Path_cnvs/$a[0]`;
            `rm $temp_dir/t.$rr $temp_dir/t.$rr.*.txt`;
        }
    }@CNVS;
}

sub smooth
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open (CHR,"$infile");
    
    # first split into chromosomes
    my $r = <CHR>;
    my $i = 0;
    my @CC;
    chomp($r);
    my @a = split("\t",$r);
    $CC[$i][0] = $a[3];
    $CC[$i][1] = $a[1]-$a[2];
    $CC[$i][2] = $a[4]+($a[5]-$a[4])/2;
    my $old_chr = $a[3];
    $i++;
    while($r = <CHR>) 
    {
        chomp($r);
        if($r =~ /NA/)
        {
            next;
        }
        @a = split("\t",$r);
        if($a[3] ne $old_chr)
        {
            dumpit(\@CC,"$outfile");
            @CC = ();
            $i = 0;
        }
        $CC[$i][0] = $a[3];
        $CC[$i][1] = $a[1] - $a[2];
        $CC[$i][2] = $a[4] + ($a[5] - $a[4]) / 2;
        $old_chr = $a[3];
        $i++;
    }
    dumpit(\@CC,"$outfile");
    close CHR;
}

##########################################smooth subs##################################
sub dumpit
{
    my @A = $_[0];
    my $outfile = $_[1];
    open (CHRO,">>$outfile");
    
    my $ref = $A[0];
    my @de_ref = @$ref; 
    my $ll = @$ref;
    my $w = 500;
    my @diff = ();
    
    # make diff array
    for(my $i = 0;$i < $ll;$i++)
    {
        push(@diff,$de_ref[$i][1]);
    }
    
    # smooth
    for(my $i = 0;$i < $ll;$i++)
    {
        my $ss = $i - 250;
        my $ee = $i + 250;
        
        if($ss < 0)
        {
            $ss = 0;
        }
        if($ee > ($ll - 1))
        {
            $ee = $ll - 1;
        }
        
        my $mm = avg(@diff[$ss..$ee]);
        
        print CHRO $de_ref[$i][0], "\t", $de_ref[$i][2], "\t", $mm, "\n";
    }
    close CHRO;
}

sub avg
{
    my(@B) = @_;
    my $ll = 0;
    my $sum = 0;
    
    for(my $i = 0;$i < scalar(@B);$i++)
    {
        if($B[$i] =~ /^[\-\.0-9]+$/)
        {
            $sum += $B[$i];
            $ll++;
        }
    }
    if ($sum == 0)
    {
        return 0;
    }
    else
    {
        return($sum / $ll);
    }
}
##########################################end of smooth subs##################################

sub delin_cnv
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my ($outfile,$cut) = @_;
    
    open (DCI,"$infile");
    open (DCO,">>$outfile");
    # first split into chromosomes
    my $r = <DCI>;
    my $in = 0;
    chomp($r);
    my @a = split("\t",$r);
    $a[1] =~ s/\.5$//;#Remove .5 of middle position and make it as integer;
    my $old_chr = $a[0];
    my $old_pos = $a[1];
    my $ss=0;
    if(abs($a[2]) >= $cut)
    {
        $in = 1;
        $ss = $old_pos;
    }
    
    while($r = <DCI>) 
    {
        chomp($r);
        @a = split("\t",$r);
        $a[1] =~ s/\.5$//;
        if(($a[0] ne $old_chr) & ($in == 1))
        {
            print DCO $old_chr, "\t", $ss, "\t", $old_pos, "\n"; $in=0; $ss=0;
        }
        
        if(($in == 1) & (abs($a[2]) < $cut))
        { # leaving zone
            print DCO $old_chr, "\t", $ss, "\t", $old_pos, "\n"; $in=0;
        }
        if((abs($a[2]) >= $cut) & ($in == 0))
        {
            $in = 1;
            $ss = $a[1];
        } # entering zone
        $old_chr = $a[0];
        $old_pos = $a[1];
    }
    close DCI;
    close DCO;
}

1;
