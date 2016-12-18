package TCGA_Lib::Bad_SNPS_and_CNVS;

use FindBin qw($Bin);
use strict;
use warnings;
use lib "$Bin/";
use MCE::Map;
use Parsing_Routines;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

my $parsing = TCGA_Lib::Parsing_Routines->new;

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
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(PB,$infile) or die("Can't open file for input: $!\n");
    open(PBO,">$outfile") or die("Can't open file for output: $!\n");
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
    close(PB);
    close(PBO);
}

sub mv_bcode
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(MVB, $infile) or die("Can't open file for input: $!");
    open(MVBO, ">$outfile") or die("Can't open file for output: $!");
    while(my $r = <MVB>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print MVBO $r, "\t", substr($a[0],0,12), "\t", substr($a[0],13,2), "\n";
    }
    close(MVB);
    close(MVBO);
}

sub pull_normal
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(PN,$infile) or die("Can't open file for input: $!\n");
    open(PNO,">$outfile") or die("Can't open file for output: $!\n");
    #removes header line;
    my $r = <PN>;
    while($r = <PN>)
    {
        chomp($r);
        my @a = split("-",$r);
        my $tumor_type = $a[-1];
        chop($tumor_type);
        my @id = split(":",$r);
        my $tid = $id[-1];
        chop($tid);
        my $nid = substr($tid,0,12);
        if($tumor_type > 9)
        {
            print PNO $r, "\t", $tid, "\t", $nid, "\n";
        }
    }
    close(PN);
    close(PNO);
}

sub hv_cd
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(HVC,$infile) or die("Can't open file for input: $!\n");
    open (HVCO,">$outfile") or die("Can't open file for output: $!\n");
    while(my $r = <HVC>)
    {
        chomp($r);
        my @a = split("\t",$r);
        unless($a[3] eq 'NaN')
        {
            print HVCO $r, "\n";
        }
    }
    close(HVC);
    close(HVCO);
}

sub dump_non_01_10_11
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(DNI,$infile) or die("Can't open file for input\n");
    open(DNO,">$outfile") or die("Can't open file for output\n");
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
    close(DNI);
    close(DNO);
}

sub mk_tn_tables
{
    my $snp_raw_files_dir = shift;
    my $self = $snp_raw_files_dir and $snp_raw_files_dir = shift if ref $snp_raw_files_dir;
    my $raw_file_name = shift;
    
    `ls "$snp_raw_files_dir" > "$raw_file_name"\_raw_files`;

    open(HC,"$raw_file_name\_raw_files") or die "Can't open file: $!\n";
    open(HCT,">$raw_file_name\_T") or die "Can't open file: $!\n";
    open(HCN,">$raw_file_name\_N") or die "Can't open file: $!\n";
    
    while(my $r = <HC>)
    {
        chomp($r);
        my ($UUID,$TCGA) = split(":",$r);
        my @tcga_array = split("-",$TCGA);
        my $id = join("-",@tcga_array[0..2]);
        my $tumor_type = $tcga_array[-1];
        $tumor_type =~ s/\D+//;
        if ($tumor_type > 9)
        {
            print HCN $TCGA, "\t", $r, "\t", $id, "\t", $tumor_type, "\n"; 
        }
        elsif($tumor_type <= 9)
        {
            print HCT $TCGA, "\t", $r, "\t", $id, "\t", $tumor_type, "\n";   
        }
        else
        {
            print "Something is wrong with $tumor_type\n";
            exit;
        }
    }
    close(HCN);
    close(HCT);
}

sub run_snps
{
    my $raw_snp_dir = shift;
    my $self = $raw_snp_dir and $raw_snp_dir = shift if ref $raw_snp_dir;
    my $infile = shift;
    my $RNA_Path = shift;
    
    my ($bird,$copy) = split(",",$raw_snp_dir);
    $bird =~ s/\/$//;
    $copy =~ s/\/$//;
    mce_map_f
    {
        chomp($_);
        my @a = split("\t",$_);
        print STDERR "The job was submitted: $a[3]\n";# ? $a[3] is TCGA ID;
        my $tmp = $parsing->temp_filename();
        if (-f "$bird/$a[1]" and -f "$bird/$a[2]")
        {
            $parsing->vlookup("$bird/$a[2]",1,"$bird/$a[1]",1,2,"y","$tmp.txt");
            flag_snps("$tmp.txt","$RNA_Path/$a[3]");
        }
        elsif(-f "$copy/$a[1]" and -f "$copy/$a[2]")
        {
            $parsing->vlookup("$copy/$a[2]",1,"$copy/$a[1]",1,2,"y","$tmp.txt");
            flag_snps("$tmp.txt","$RNA_Path/$a[3]"); 
        }
        else
        {
            print "SNP array data not available for $a[1].\n";
        }
        `rm -f $tmp.txt`;
    }$infile; 
}

sub flag_snps
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(FSI,$infile) or die("Can't open file for input: $!\n");
    open(FSO,">$outfile") or die("Can't open file for output: $!\n");
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
    close(FSI);
    close(FSO);
}

sub get_cd_bed
{
    my $affy_cd = shift;
    my $self = $affy_cd and $affy_cd = shift if ref $affy_cd;
    my $RNA_Path = shift;
    my $bad_snps = shift;
    
    opendir(DD,"$bad_snps") or die "no $bad_snps: $!\n";
    my @dd = readdir(DD);
    closedir(DD);
    @dd = grep{!/\.$/}@dd;#gets rid off . and ..
    mce_map
    {
      print STDERR "working on $_\n";
        
      my $tmp = $parsing->temp_filename();
       `sort -u -k 1,1 $bad_snps/$_ > $bad_snps/$_.uniq`;
       $parsing->vlookup("$bad_snps/$_.uniq",1,"$affy_cd",4,"1,2,3","y","$bad_snps/$tmp.1.txt");
       `grep NaN -v $bad_snps/$tmp.1.txt > $bad_snps/$tmp.2.txt`;
       $parsing->pull_column("$bad_snps/$tmp.2.txt","2,3,4","$bad_snps/$tmp.3.txt");
       `sort -k 1,1 -k 2,2n $bad_snps/$tmp.3.txt > $RNA_Path/bad_snps_bed/$_.bed`;
       `rm $bad_snps/$tmp.*.txt`;
    }@dd;
}

sub run_cnv
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $copy = shift;#directory where SNP_Array data for copynumber data is
    my $bird = shift;#directory where SNP-Array data for birdiseed is
    my $cnv_bed = shift;#./affy6/cnv.hg19.bed
    my $RNA_Path_cnvs = shift;
    
    `mkdir temp` unless(-d "temp");

    open(my $CIN,"$infile") or die "Can't open $infile: $!\n";
    my @CNVS = <$CIN>;

    mce_map
    {
        chomp($_);
        my @a = split("\t",$_);
        print STDERR "The job was submitted: $a[3]\n";
        
        my $rr = rand();
        $rr = substr($rr,2,6);
        #This will check if the CNV files exist. They should always exist unless there were changes made to the directory prior to this running.
        if (-e "$copy/$a[5]" and -e "$copy/$a[6]")
        {
            $parsing->vlookup("$copy/$a[5]",1,"$copy/$a[6]",1,2,"y","temp/t.$rr");
            $parsing->vlookup("temp/t.$rr",1,$cnv_bed,4,"1,2,3","y","temp/t.$rr.1.txt");
            $parsing->strip_head("temp/t.$rr.1.txt","temp/t.$rr.2.txt",2);
            `sort -k 4,4 -k 5,5n temp/t.$rr.2.txt > temp/t.$rr.3.txt`;
            smooth("temp/t.$rr.3.txt","temp/t.$rr.4.txt");
            delin_cnv("temp/t.$rr.4.txt","temp/t.$rr.5.txt",0.5);
            `grep NaN -v temp/t.$rr.5.txt > $RNA_Path_cnvs/$a[3].bed`;
            `rm temp/t.$rr temp/t.$rr.*.txt`;
        }
    }@CNVS;
}

sub smooth
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(CHR,"$infile") or die "Can't open $infile: $!\n";
    
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
    open(CHRO,">>$outfile") or die "Can't open $outfile: $!\n";
    
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
########################################################################

sub delin_cnv
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    my $cut = shift;
    
    open(DCI,"$infile") or die "Can't open $infile: $!\n";
    open(DCO,">>$outfile") or die "Can't open $outfile: $!\n";
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