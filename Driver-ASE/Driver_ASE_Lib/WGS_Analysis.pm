package Driver_ASE_Lib::WGS_Analysis;

use strict;
use warnings;
use FileHandle;
use FindBin qw($Bin);
use lib "$Bin/";
use Parsing_Routines;
use autodie;
no warnings 'once';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

my $parsing = Driver_ASE_Lib::Parsing_Routines->new;

###############################################################################
# WGS_Analysis.pm version 1
# Written by Tomas Babak, Zhongshan Cheng, Michael Vermeulen and Micheal Rollins-Green
#  in the Department of Biology at Queen's University.
#
# WGS_Analysis
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

###############################subs#############################\

sub filter_not_done_somatic
{
    my $self = shift;
    my ($not_done,$lookup_file,$compare_file,$WGS_Path,$somatic) = @_;
    
    $parsing->vlookup("$lookup_file",2,"$compare_file",1,1,"y","$WGS_Path/$somatic\_look.txt");
    my @filter_done = `cat $WGS_Path/$somatic\_look.txt`;
    
    open (SO,">$not_done");
    for my $sl (@filter_done)
    {
        if ($sl =~ /\tNaN+$/i)
        {
            chomp($sl);
            my $mp = [split("\t",$sl)]->[0];
            print SO $mp,"\n";
        }
    }
    close (SO);
}

sub Varscan_filter
{
    my $not_done_wgs_som_file = shift; #path to the directory where wgs mpileups are stored
    my $self = $not_done_wgs_som_file and $not_done_wgs_som_file = shift if ref $not_done_wgs_som_file;
    my $somatic = shift; #directory where somatic data will be stored
    my $readcutoff = shift || 20;
    my $VarType = shift || ".*"; #what will be filtered in this routine (Somatic)
    my $normal_alt_frq = shift || 1;
    $normal_alt_frq *= 100;
    my $tumor_alt_frq = shift || 0;
    $tumor_alt_frq *= 100;

    #Columns of Varscan output;
    #ARRAY_NUM: HEADERS           ROWELMENTS
    #0: chrom                     1
    #1: position                  10250
    #2: ref                       A
    #3: var                       C
    #4: normal_reads1             9
    #5: normal_reads2             0
    #6: normal_var_freq           0%
    #7: normal_gt                 A
    #8: tumor_reads1              7
    #9: tumor_reads2              3
    #10: tumor_var_freq           30%
    #11: tumor_gt                 M
    #12: somatic_status           Somatic
    #13: variant_p_value          1.0
    #14: somatic_p_value          0.1238390092879262
    #15: tumor_reads1_plus        2
    #16: tumor_reads1_minus       5
    #17: tumor_reads2_plus        1
    #18: tumor_reads2_minus       2
    #19: normal_reads1_plus       3
    #20: normal_reads1_minus      6
    #21: normal_reads2_plus       0
    #22: normal_reads2_minus      0
    my @dd;
    #directory or file will be passed to routine based on if user wants to overlap data or not

    @dd = `cat $not_done_wgs_som_file`;
    
    for (my $i = 0;$i < scalar(@dd);$i++)
    {
        chomp ($dd[$i]);
        my $newfile = $dd[$i];
        if ($newfile =~ /\//)
        {
            $newfile = [split("/",$newfile)]->[-1];
        }
        $newfile =~ s/\.(snp|indel)//;
        $newfile =~ s/\.varscan//;
        
        open (my $V,"$dd[$i]");
        
        open (SOM,">>$somatic/$newfile");
        while (my $l = <$V>)
        {
            chomp($l);
            next if $l =~ /^chrom/;
            next if $l =~ /^(GL|MT)/;
            next if $l =~ /^\s*$/;
            my @as = split("\t",$l);
            if (!defined $as[5])
            {
                next;
            }
            my $N_n = $as[4] + $as[5];
            my $T_n = $as[8] + $as[9];
            my $T_alt_frq = $as[10];           
            $T_alt_frq =~ s/%//g;
            
            my $N_alt_frq = $as[6];
            $N_alt_frq =~ s/%//g;
            #only consider specific type variants;   
            next unless ($as[12] =~ /$VarType/i);
            my $print_value;
            if ($N_n >= $readcutoff and $T_n >= $readcutoff and $N_alt_frq <= $normal_alt_frq and $T_alt_frq >= $tumor_alt_frq and $as[3] !~ /[+-]/)
            {
                #for single somatic mutations
                $print_value = VarscanTable2Simple($l);
                print SOM $print_value,"\n" if $l =~ /^(chr)*(\d+|x|y)\b/i;
                
            }
            elsif ($as[3] =~ /[+-]/ and $N_n >= $readcutoff and $T_n >= $readcutoff)
            {
                #for indels, only consider read cutoff;
                $print_value = VarscanTable2Simple($l);
                print SOM $print_value,"\n" if $l =~ /^(chr)*(\d+|x|y)\b/i;
            }
            
            #For het geno in normal, that is T_alt_frq >45;
            #By manually checking, only about 0.1% percent are in this case;
            #and most of somatic muations called by varscan are actually SNPs;
            #Thus no need to consider further!
        }
        close ($V);
        close (SOM);
    } 
}

sub VarscanTable2Simple
{
    my $l = shift;
    my @as = split("\t",$l);
    $as[0] = "chr".$as[0] unless $as[0] =~ /^chr/i;
    my $ref_alt = $as[2]."|".$as[3];
    my $out = $as[0]."-".$as[1]."-".$ref_alt."\t".$as[0]."\t".($as[1]-1)."\t".$as[1]."\t".$ref_alt."\t1";
    return $out;
}

sub mk_files_for_wgs
{
    my $somatic = shift; 
    my $self = $somatic and $somatic = shift if ref $somatic;
    my $somatic_list = shift;
    
    my @dd;
    if (-d $somatic)
    {
        opendir (DD,$somatic);
        @dd = grep{!/^\./ && -f "$somatic/$_"} readdir(DD);
        closedir (DD);
    }
    elsif (-f $somatic)
    {
        @dd = `cat $somatic`;
    }
    
    open (MFO,">$somatic_list");
   
    for (my $i = 0;$i < scalar(@dd);$i++)
    {
        chomp($dd[$i]);
        if (-d $somatic)
        {
            print MFO $somatic."/".$dd[$i], "\t", $dd[$i], "\n";
        }
        elsif (-f $somatic)
        {
            my $somatic_file = [split("/",$dd[$i])]->[-1];
            print MFO $dd[$i],"\t",$somatic_file,"\n";
        }
    }
    close (MFO);
}

sub var_ID_to_bed
{
    my $rowlabels = shift; #file that contains rowlabels for wgs
    my $self = $rowlabels and $rowlabels = shift if ref $rowlabels;
    my $varsbed = shift; #output file for rowlabels
    
    open (VTBI,"$rowlabels");
    open (VTBO,">$varsbed");
    while (my $r = <VTBI>)
    {
        chomp($r);
        my @a = split("\-",$r);
        
        print VTBO $a[0], "\t", $a[1]-1, "\t", $a[1], "\t", $r, "\n";
    }
    close (VTBI);
    close (VTBO);
}

sub simpl
{
    my $reg_file = shift; #path to files from reg directory in the Database directory
    my $self = $reg_file and $reg_file = shift if ref $reg_file;
    my $output_bed = shift; #file to output data
    
    open (SI,"$reg_file");
    open (SO,">$output_bed");
    while (my $r = <SI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print SO join("\t",@a[0..2]), "\t", 1 , "\n";
    }
    close (SI);
    close (SO);
}

sub up_down_tss
{
    my $bed_file = shift; #input bed file
    my $self = $bed_file and $bed_file = shift if ref $bed_file;
    my ($upstream_of_tss_len,$downstream_of_tss_len,$output_bed) = @_; #upstream number, downstream number, bed file to output results
    
    open (UDTI,"$bed_file");
    open (UDTO,">$output_bed");
    while (my $r = <UDTI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        
        if ($a[5] eq '+')
        {
            if ($a[1] > $upstream_of_tss_len)
            {
                print UDTO $a[0], "\t", $a[1] - $upstream_of_tss_len, "\t", $a[1] + $downstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
            }
            else
            {
                print  UDTO $a[0], "\t", 0, "\t", $a[1] + $downstream_of_tss_len,  "\t", join("\t",@a[3..5]),"\n";
            }
        }
        
        if ($a[5] eq '-')
        {
            if ($a[2] > $downstream_of_tss_len)
            {
                print UDTO $a[0], "\t", $a[2] - $downstream_of_tss_len, "\t", $a[2] + $upstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
            }
            else
            {
                print  UDTO $a[0], "\t", 0, "\t", $a[2] + $upstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
            }
        }
    }
    close (UDTI);
    close (UDTO);
}

sub print_1
{
    my $upstream_bed = shift; #bed file from up_down_tss
    my $self = $upstream_bed and $upstream_bed = shift if ref $upstream_bed;
    my $output_bed = shift; #bed file that will have data written to it
    
    open (PRI,"$upstream_bed");
    open (PRO,">$output_bed");
    while (my $r = <PRI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print PRO join("\t",@a[0..2]), "\t", 1, "\n";
    }
    close (PRI);
    close (PRO);
}

sub up_tss_gene
{
    my $upstream_bed = shift; #bed file from up_down_tss
    my $self = $upstream_bed and $upstream_bed = shift if ref $upstream_bed;
    my ($upstream_of_tss_len,$output_bed) = @_; #upstream length number, bed file that will have data written to it
    
    open (UTGI,"$upstream_bed");
    open (UTGO,">$output_bed");
    while (my $r = <UTGI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        
        if ($a[5] eq '+')
        {
            if ($a[1] > $upstream_of_tss_len)
            {
                print UTGO $a[0], "\t", $a[1] - $upstream_of_tss_len, "\t", $a[2], "\t", join("\t",@a[3..5]), "\n";
            }
            else
            {
                print UTGO $a[0], "\t", 0, "\t", $a[2], "\t", join("\t",@a[3..5]),"\n";
            }
        }
        
        if ($a[5] eq '-')
        {
            print UTGO $a[0], "\t", $a[1], "\t", $a[2] + $upstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
        }
    }
    close (UTGI);
    close (UTGO);
}

sub vlookem_all
{
    my $overlap = shift; #path to the directory that holds beds made from routine print_1
    my $self = $overlap and $overlap = shift if ref $overlap;
    my ($vars_bed,$an_vars,$overlapSelect,$stream_data) = @_; #path to output bed file, path to directory where annotation data is stored, overlapSelect command, name of stream data file
    
    opendir (DD,$overlap);
    my @dd = grep{!/^\./ && -f "$overlap/$_"} readdir(DD);
    closedir (DD);
    
    my $cc = 1;
    
    open (COL,">collabels.txt");
    
    for (my $i = 0;$i < scalar(@dd);$i++)
    {
        my $cc2 = $cc + 1;
        
        print "Working on $dd[$i]\n";
        print COL $dd[$i], "\n";
        
        `$overlapSelect $overlap/$dd[$i] $vars_bed -idOutput out`; 
        $parsing->vlookup("$an_vars/$stream_data"."$cc",1,"out",1,2,"y","$an_vars/$stream_data"."$cc2");
        $cc++;
    }
    close (COL);
}

sub my_up_down_tss
{
    my $refseq_hg9_bed = shift; #path to bed refseq bed file in Database directory
    my $self = $refseq_hg9_bed and $refseq_hg9_bed = shift if ref $refseq_hg9_bed;
    my ($upstream_of_tss_len,$downstream_of_tss_len,$refseq_hg9_output) = @_; #upstream length number, downstream length number, file to output results
    
    open (MUDTI,"$refseq_hg9_bed");
    open (MUDTO,">$refseq_hg9_output");
    while (my $r = <MUDTI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        
        if ($a[5] eq '+')
        {
            if ($a[1] > $upstream_of_tss_len)
            {
                #be cautious of $upstream of tss length
                print MUDTO $a[0], "\t", $a[1]-$upstream_of_tss_len, "\t", $a[2], "\t", join("\t",@a[3..5]), "\n";
            }
            else
            {
                print MUDTO $a[0], "\t", 0, "\t", $a[2],  "\t", join("\t",@a[3..5]),"\n";
            }
        }
        
        if ($a[5] eq '-')
        {
            if ($a[2] > $downstream_of_tss_len)
            {
                #be cautious of $downstream of tss length
                print MUDTO $a[0], "\t", $a[1], "\t", $a[2] + $upstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
            }
            else
            {
                print  MUDTO $a[0], "\t", 0, "\t", $a[2] + $upstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
            }
        }
    }
    close (MUDTI);
    close (MUDTO);
}

sub mk_coding
{
    my $var_bed = shift; #bed file create from var_ID_to_bed routine
    my $self = $var_bed and $var_bed = shift if ref $var_bed;
    my $var_bed_output = shift; #file to output results
    
    open (MCI,"$var_bed");
    open (MCO,">$var_bed_output");
    
    while (my $r = <MCI>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        
        $a[3] =~ /^([^\-]+)\-([^\-]+)\-(.*)$/;
        
        my $tt = $3;
        
        $tt =~ s/\|/,/;
        my @w2 = split(",",$tt);
        
        if ($w2[1] =~ /[+\-]/)
        {
            
            my $t = $w2[1];
            
            $t =~ s/^[\+\-]//;
            
            my $ll = length($t)/3;
            if ($ll =~ /\./) # no triplet - disruptive
            {
               $a[3] = $w2[0]."|"."frameshift";        
            }
            else
            {
               $a[3] = $w2[0]."|"."triplet_shift";
            }
        }
        else
            {
                $a[3] = join("\|",@w2);
            } 
        print MCO join("\t",@a), "\n";
    }
    close (MCI);
    close (MCO);
}

sub Syn
{
    my $subs_bed = shift; #input bed file with default name subs.bed
    my $self = $subs_bed and $subs_bed = shift if ref $subs_bed;
    my ($genome,$subs_syn,$dbpath,$overlapSelect,$ccds_bed,$ccds_tab) = @_; #genome hg19, file to output data with default name subs.syn, path to the Database directory, overlapSelect command, bed file in the Database directory with the default name of ccds.hg19.bed, tab file in the Database directory with the default name of ccds.hg19.tab
    
    open (SO,">$subs_syn");
    
    my $syn = 0;
    my $nsyn = 0;
    
    #data format of ccds
    # 0: chr1
    # 1: 67000041
    # 2: 67208778
    # 3: CCDS30744.1
    # 4: 0
    # 5: +
    # 6: 67000041
    # 7: 67208778
    # 8: 0
    # 9: 25
    # 10: 10,64,25,72,57,55,176,12,12,25,52,86,93,75,501,128,127,60,112,156,133,203,65,165,23,
    # 11: 0,91488,98711,101585,105418,108451,109185,126154,133171,136636,137585,138922,142645,145319,147510,154789,155831,161075,184935,194905,199389,204976,206299,206913,208714,
    
    if ($genome eq 'hg19')
    {
        $ccds_bed = "$dbpath/$ccds_bed";
        print STDERR "$ccds_bed does not exist in directory $dbpath!\n" and exit unless(-e $ccds_bed);
        $ccds_tab = "$dbpath/$ccds_tab";
        print STDERR "$ccds_tab does not exit in directory $dbpath!\n" and exit unless(-e $ccds_tab);    
    }
    else
    {
        print STDERR "Genome is either unknown or is not used within this pipeline.\n";
        exit;
    }

    print "Overlapping with $ccds_bed\n";
    my $select_command = "$overlapSelect $ccds_bed $subs_bed ccds.input.out -mergeOutput";
    print $select_command,"\n";
    `$select_command`;
    
    #data format of ccds.input.out;
    # 0: chr1
    # 1: 1431164
    # 2: 1431165
    # 3: C|T
    # 4: chr1
    # 5: 1407264
    # 6: 1431197
    # 7: CCDS30.1
    # 8: 0
    # 9: +
    # 10: 1407264
    # 11: 1431197
    # 12: 0
    # 13: 16
    # 14: 205,77,102,60,70,166,70,156,57,126,125,52,71,168,109,333,
    # 15: 0,5389,6759,7164,8982,10253,10660,13131,13897,14225,14659,15978,17319,18372,18678,23600,
    
    open (F,"ccds.input.out");
    
    #make lookup
    print "Making sequence lookup table\n";
    my %hash = ();
    open (LOOK,$ccds_tab);
    while (my $r = <LOOK>)
    {
        $r =~ /(.*?)\t(.*?)\n/;#hash of ccds id and fasta seq;
        $hash{$1} = $2;
    }
    
    print "Going through your file\n";
    while (my $r = <F>)
    {
        my @a = split("\t",$r);
        
        # what are the two observed bases?
        my $aa;
        my $bb;
        if ( $a[3] =~ /([ACGT])\|([ACGT])/)
        {
            $aa = $1;
            $bb = $2;
        }
        else
        {
            print $a[3], " does not have any base identities\n";
            next;
        }
        
        if ($aa eq '0')
        {
            next;
        }
        if ($a[9] eq '-')
        { 
            $aa =~ tr/ACGT/TGCA/; 
            $bb =~ tr/ACGT/TGCA/;
        }
        # print Dominant base and ratio of canonical:non-canonical counts 
        
        chomp($a[15]);
        # what is my codon?
        
        # Return three nt of amino acid harboring mutation
        # mutation pos in the three nt of amino acid
        # and the codon start pos
        #$a[1]:allele start pos;
        #$a[5]:gene start pos;
        #$a[6]:gene end pos;
        #$a[9]:strand;
        #$a[14]:exon length;
        #$a[15]:exon start;
        if (exists $hash{$a[7]})
        {
            my($codon,$pos,$cc_start) = get_codon($a[1],$a[5],$a[6],$a[9],$a[14],$a[15],$hash{$a[7]});
            
            #what kind of substitution is it?
            my $can_aa = codon($codon);
            my $new_aa = codon(substr($codon,0,$pos).$bb.substr($codon,$pos+1,3-$pos-1));
            
            # severity 0=non-synonymous, 1=missense(same group) 2=missense(different group) 3=stop codon(nonsense)
            my $severity = severity($can_aa,$new_aa);
            
            print SO $a[0], "\t", $a[1], "\t", $a[2], "\t", $a[3],
            "\t", $a[7], "\t", $a[8], "\t", $a[9], "\t", $codon, "\t",
            substr($codon,0,$pos).$bb.substr($codon,$pos+1,3-$pos-1),
            "\t", $can_aa, "\t", $new_aa, "\t", $severity,
            "\t",($cc_start/3+1), "\n";
            
            if ($can_aa eq $new_aa)
            {
                $syn++;
            }
            else
            {
                $nsyn++;
            }
        }
        else
        {
            next;
        }
    }
    
    print "Total number of events: ", $nsyn + $syn, "\n";
    print "Synonymous substitutions: $syn\n";
    print "NonSynonymous substitutions: $nsyn\n";
    if ($syn == 0 && $nsyn == 0)
    {
        print "Synonymous/NonSynonymous: 0\n";
    }
    else
    {
       print "Synonymous/NonSynonymous: ", sprintf "%0.*f", 2,$syn/$nsyn;
       print "\n";
    }
    
    `rm ccds.input.out`;
    close (SO);
}

#########################################end of Syn subs#######################################
sub get_codon
{
    my($poss,$s,$e,$orient,$lens,$starts,$seq) = @_;

    # get start within transcript-go through one exon at a time, when within, how far to start?
    # which exon is it in?
    $poss -= $s;
    
    my @ll = split(",",$lens);#exon length;
    my @ss = split(",",$starts);#exon start pos;
    push(@ss,$e - $s);
    my $txs = 0;
    
    for (my $i = 1;$i < scalar(@ss);$i++)
    {
        if (($poss >= $ss[$i-1]) & ($poss < $ss[$i]))
        {
            #add up exon lenghts to here
            for (my $j = 0;$j < $i-1;$j++)
            {
                $txs += $ll[$j];
            }
            #add up the exon length from the start of the next exon to muation pos;	   
            $txs += ($poss - $ss[$i-1]);
            last;#stop further searching;    
        } # if in this exons
    } #end for
    
    # if orient backwards, invert start
    if ($orient eq '-')
    {
        my $txlen = 0;
        for (my $i = 0;$i < scalar(@ll);$i++)
        {
            $txlen += $ll[$i];#Get the total length of exons;
        } #end for
        $txs = $txlen - $txs - 1;#Get the length from exon start to mutation pos;
    }
    
    # now get the right codon and spit out the position
    # based on that fact that one amino acid is made up by THREE nt;
    # get the codon start nt and mutation pos in the three nt of amino acid
    my($codon_start,$codon_pos) = getPosFromPos($txs);
    
    # Return three nt of amino acid
    # mutation pos in the three nt of amino acid
    # and the codon start pos
       return(substr($seq,$codon_start,3),$codon_pos,$codon_start); 
    
}

sub getPosFromPos
{
    my($s) = @_;
    
    my $x = $s/3;
    
    $x =~ /([0-9]+)\.*(.*)/;
    my $start = $1*3;
    my $duss = $2;
    if ($duss eq '')
    {
        return($start,0);
    }#perfect matching the end of ABC;
    elsif ($duss =~ /3/)
    {
        return($start,1);
    }#matching the nt after ABCa;
    else
    {
        return($start,2);
    }#matching the second nt after ABCab;
}

sub codon
{

    my($a) = @_;
    
    my $codontable = 'ATT I
    ATC	I
    ATA	I
    CTT	L
    CTC	L
    CTA	L
    CTG	L
    TTA	L
    TTG	L
    GTT	V
    GTC	V
    GTA	V
    GTG	V
    TTT	F
    TTC	F
    ATG	M
    TGT	C
    TGC	C
    GCT	A
    GCC	A
    GCA	A
    GCG	A
    GGT	G
    GGC	G
    GGA	G
    GGG	G
    CCT	P
    CCC	P
    CCA	P
    CCG	P
    ACT	T
    ACC	T
    ACA	T
    ACG	T
    TCT	S
    TCC	S
    TCA	S
    TCG	S
    AGT	S
    AGC	S
    TAT	Y
    TAC	Y
    TGG	W
    CAA	Q
    CAG	Q
    AAT	N
    AAC	N
    CAT	H
    CAC	H
    GAA	E
    GAG	E
    GAT	D
    GAC	D
    AAA	K
    AAG	K
    CGT	R
    CGC	R
    CGA	R
    CGG	R
    AGA	R
    AGG	R
    TAA	Stop
    TAG	Stop
    TGA	Stop
    ';
    
    $codontable =~ /$a\t(.*)/;
    
    return $1;
}

sub severity 
{
    my($can_aa,$new_aa) = @_;
    
    if ($can_aa eq $new_aa)
    {
        return 0;
    }
    if ($new_aa eq 'Stop')
    {
        return 3;
    }
    
    my $hydrophobic = 'AILMFWYV';
    my $polar = 'STNQ';
    my $positive = 'RHK';
    my $negative = 'DE';
    
    # the rest of them are weird so a sub to these = missense 2
    
    if (($hydrophobic =~ /$can_aa/)&($hydrophobic =~ /$new_aa/))
    {
        return 1;
    }
    if (($polar =~ /$can_aa/)&($polar =~ /$new_aa/))
    {
        return 1;
    }
    if (($positive =~ /$can_aa/)&($positive =~ /$new_aa/))
    {
        return 1;
    }
    if (($negative =~ /$can_aa/)&($negative =~ /$new_aa/))
    {
        return 1;
    }
    
    return 2;
}
############################End of Syn subs#############################################

sub mk_uniq_label
{
    my $subs_syn = shift; #syn file made from routine Syn
    my $self = $subs_syn and $subs_syn = shift if ref $subs_syn;
    my $subs_syn_out = shift; #file that will be written to
    
    open (MULI,"$subs_syn");
    open (MULO,">$subs_syn_out");
    while (my $r = <MULI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print MULO $a[0]."-".$a[1], "\t", $r, "\n";
    }
    close (MULI);
    close (MULO);
}

sub process_syn
{
    my $subs_syn_process = shift; #syn file that will be processed
    my $self = $subs_syn_process and $subs_syn_process = shift if ref $subs_syn_process;
    my $subs_process_out = shift; #file that will have written results
    
    open (PSI,"$subs_syn_process");
    open (PSO,">$subs_process_out");
    while (my $r = <PSI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print PSO $a[0]."-".$a[1], "\t", $a[9].$a[12].$a[10],
        "\t", $a[11], "\n";
    }
    close (PSI);
    close (PSO);
}

sub process_shifts
{
    my $overlap_input = shift; #file with results from overlapSelect command
    my $self = $overlap_input and $overlap_input = shift if ref $overlap_input;
    my $shifts_append = shift; #file that will have results appended to
    
    open (PRSI,"$overlap_input");
    open (PRSO,">>$shifts_append");
    while (my $r = <PRSI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        $a[3] =~ s/\|/,/;
        my @w=split(",",$a[3]);
        my $t = $w[1];
        if ($t eq 'frameshift')
        {
            print PRSO $a[0]."-".$a[1], "\t", 'frameshift', "\t", 3, "\n";
        }
        elsif ($t eq 'triplet_shift')
        {
            print PRSO $a[0]."-".$a[1], "\t", 'triplet_shift', "\t", 2, "\n";
        }
        else
        {
            print "shleisse \n";
            next;
        }
    }
    close (PRSI);
    close (PRSO);
}

sub pt
{
    my $rowlabels = shift; #rowlabels file from matricize routine in the Parsing_Routines module
    my $self = $rowlabels and $rowlabels = shift if ref $rowlabels;
    my $rowlabels_out = shift; #file that will have results written to it
    
    open (PTI,"$rowlabels");
    open (PTO,">$rowlabels_out");
    while (my $r = <PTI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        $a[0] =~ /^([^\-]+)\-([^\-]+)\-/;
        
        my $xx=$2-1;
        
        print PTO $1."-".$xx, "\t", $r, "\n";
    }
    close (PTI);
    close (PTO);
}

1;
