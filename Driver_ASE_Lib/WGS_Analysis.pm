package Driver_ASE_Lib::WGS_Analysis;

use strict;
use warnings;
use FileHandle;
use FindBin qw($Bin);
use lib "$Bin/";
use Parsing_Routines;

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

###############################subs#############################
sub launch_prune
{
    my $indir = shift;#mpileups_directory
    my $self = $indir and $indir = shift if ref $indir;
    my $outdir = shift;#output_directory
    my $MAN = shift;#max_alt_n
    my $MAT = shift;#min_alt_t
    my $MC = shift;#min_cov
    
    opendir(DD,$indir) or die "no dir\n";
    my @dd = readdir(DD);
    closedir(DD);
    
    for(my $i = 0;$i < scalar(@dd);$i++)
    {
        if($dd[$i] =~ /^\.+$/)
        {
            next;
        }
        $dd[$i] =~ /^(.*?)\./;
        my $out = $1;
        print STDERR "processing $dd[$i]\n";
        prune_varscan("$indir/$dd[$i]",$MAN,$MAT,$MC,"$outdir/$out");
    }
}

sub prune_varscan
{
    my $indir = shift;
    my $self = $indir and $indir = shift if ref $indir;
    my $MAN = shift;
    my $MAT = shift;
    my $MC = shift;
    my $outdir = shift;
    
    open(PVI,"$indir") or die("Can't open file for output: $!\n");
    open(PVO,">$outdir") or die("Can't open for output: $!\n");
    my $r = <PVI>;

    while($r = <PVI>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        
        # only somatic calls
        unless($a[12] =~ /Somatic/gi)
        {
            next;
        }
        
        # min cov
        my $n_cov = $a[4] + $a[5];
        my $t_cov = $a[8] + $a[9];
        unless(($n_cov >= $MC) & ($t_cov >= $MC))
        {
            next;
        }
        
        # min prop alt
        unless(($a[5] / $n_cov) <= $MAN)
        {
            next;
        }
        unless(($a[9] / $t_cov) >= $MAT)
        {
            next;
        }
        
        if($a[0] =~ /^[0-9xy]+$/i)
        {
            $a[0] = "chr".$a[0];
        }
        
        my $var = $a[2]."|".alt_it($a[2],$a[11]);
        
        print PVO $a[0]."-".$a[1]."-".$var, "\t", $a[0], "\t", $a[1]-1, "\t", $a[1], "\t", $var, "\t", 1, "\n";
    }
    close(PVI);
    close(PVO);
}

sub Varscan_filter
{
    my $input = shift;
    my $self = $input and $input = shift if ref $input;
    die unless (-e "$input");
    my $output = shift;
    my $readcutoff = shift || 20;
    print STDERR "You need to provide readcutoff in numeric\n" and die if $readcutoff =~ /^\D+/i;
    my $VarType = shift || ".*";
    my $normal_alt_frq = shift || 1;
       $normal_alt_frq *= 100;
    my $tumor_alt_frq = shift || 0;
       $tumor_alt_frq *= 100;
    print STDERR "Please be noted that your parameters for Varscan_Filter are listed below:\n",
                 "VarscanOutput: $input\n","readcutoff: $readcutoff\n",
                 "VarType: $VarType\n","normal_alt_frq: $normal_alt_frq\n",
                 "tumor_alt_frq: $tumor_alt_frq\n";
                 
                 
    #Columns of Varscan output;
    #ARRAY_NUM: HEADERS           ROWEELMENTS
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
    
opendir(DD,$input) or die "no dir\n";
    my @dd = readdir(DD);
       @dd=grep{!/^\./}@dd;
    closedir(DD);
    
    for(my $i = 0;$i < scalar(@dd);$i++)
    {
        my $newfile = $dd[$i];
           chomp($newfile);
           $newfile =~ s/\.(snp|indel)//;
           $newfile =~ s/\.varscan//;
        open (my $V,"$input/$dd[$i]") or die "Can not open the file $input/$dd[$i]: $!";
        open(SOM,">>$output/$newfile") or die "Can't open $output/$newfile: $!\n";
        while (my $l=<$V>)
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
            
            if ($N_n >= $readcutoff and
                $T_n >= $readcutoff and
                $N_alt_frq <= $normal_alt_frq and
                $T_alt_frq >= $tumor_alt_frq and
                $as[3] !~ /[+-]/
                )
            {
                #for single somatic mutations
                print SOM VarscanTable2Simple($l),"\n" if $l =~ /^(chr)*(\d+|x|y)\b/i;  
                
            }
            elsif($as[3] =~ /[+-]/ and
                $N_n >= $readcutoff and
                $T_n >= $readcutoff
                )
            {
                #for indels, only consider read cutoff;
                print SOM VarscanTable2Simple($l),"\n" if $l =~ /^(chr)*(\d+|x|y)\b/i;
            }
            
            #For het geno in normal, that is T_alt_frq >45;
            #By manually checking, only about 0.1% percent are in this case;
            #and most of somatic muations called by varscan are actually SNPs;
            #Thus no need to consider further!
        }
        close $V;
        close SOM;
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
#################################prune_varscan sub##############################
sub alt_it
{
    my($aa,$bb) = @_;
    
    my $iu = 'R	AG
    Y	CT
    S	GC
    W	AT
    K	GT
    M	AC';
    
    my @I = split("\n",$iu);
    
    for(my $i = 0;$i < scalar(@I);$i++)
    {
        my @w = split("\t",$I[$i]);
        if($w[0] eq $bb)
        {# this is it!
            if(substr($w[1],0,1) eq $aa)
            {
                return(substr($w[1],1,1));
            }
            elsif(substr($w[1],1,1) eq $aa)
            {
                return(substr($w[1],0,1));
            }
            else
            {
                return('problem')
            }
        }
    }
    #This may be for parsing indel;
    if($bb =~ /^[ACGT]$/)
     {
        return($bb);
    }
    elsif($bb =~ /^\*\/([\+\-].*)/)
    {
        return($1);
    }
    else
    {
        return('problem2');
    }
}
#################################end of prune_varscan sub##############################

sub mk_files_for_wgs
{
    my $indir = shift;
    my $self = $indir and $indir = shift if ref $indir;
    my $outfile = shift;
    
    opendir(DD,$indir) or die "Can't open $indir: $!\n";
    my @dd = readdir(DD);
    @dd = grep{!/\.$/}@dd;
    open(MFO,">$outfile") or die "Can't open $outfile: $!\n";
   
    for(my $i = 0;$i < scalar(@dd);$i++)
    {  
        print MFO $indir."/".$dd[$i], "\t", $dd[$i], "\n";
    }
    close(MFO);
}

sub var_ID_to_bed
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(VTBI,"$infile") or die "Can't open $infile: $!\n";
    open(VTBO,">$outfile") or die "Can't open $outfile: $!\n";
    while(my $r = <VTBI>)
    {
        chomp($r);
        my @a = split("\-",$r);
        
        print VTBO $a[0], "\t", $a[1]-1, "\t", $a[1], "\t", $r, "\n";
    }
    close(VTBI);
    close(VTBO);
}

sub simpl
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(SI,"$infile") or die "Can't open $infile: $!\n";
    open(SO,">$outfile") or die "Can't open $outfile: $!";
    while(my $r = <SI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print SO join("\t",@a[0..2]), "\t", 1 , "\n";
    }
    close(SI);
    close(SO);
}

sub up_down_tss
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $upstream_of_tss_len = shift;
    my $downstream_of_tss_len = shift;
    my $outfile = shift;
    
    open(UDTI,"$infile") or die "Can't open $infile: $!\n";
    open(UDTO,">$outfile") or die "Can't open $outfile: $!\n";
    while(my $r = <UDTI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        
        if( $a[5] eq '+')
        {
            if($a[1] > $upstream_of_tss_len)
            {
                print UDTO $a[0], "\t", $a[1] - $upstream_of_tss_len, "\t", $a[1] + $downstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
            }
            else
            {
                print  UDTO $a[0], "\t", 0, "\t", $a[1] + $downstream_of_tss_len,  "\t", join("\t",@a[3..5]),"\n";
            }
        }
        
        if( $a[5] eq '-')
        {
            if($a[2] > $downstream_of_tss_len)
            {
                print UDTO $a[0], "\t", $a[2] - $downstream_of_tss_len, "\t", $a[2]+$upstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
            }
            else
            {
                print  UDTO $a[0], "\t", 0, "\t", $a[2] + $upstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
            }
        }
    }
    close(UDTI);
    close(UDTO);
}

sub print_1
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(PRI,"$infile") or die "Can't open $infile: $!\n";
    open(PRO,">$outfile") or die "Can't open $outfile: $!\n";
    while(my $r = <PRI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print PRO join("\t",@a[0..2]), "\t", 1, "\n";
    }
    close(PRI);
    close(PRO);
}

sub up_tss_gene
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $upstream_of_tss_len = shift;
    my $outfile = shift;
    
    open(UTGI,"$infile") or die "Can't open $infile: $!\n";
    open(UTGO,">$outfile") or die "Can't open $outfile: $!\n";
    while(my $r = <UTGI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        
        if( $a[5] eq '+')
        {
            if($a[1] > $upstream_of_tss_len)
            {
                print UTGO $a[0], "\t", $a[1] - $upstream_of_tss_len, "\t", $a[2], "\t", join("\t",@a[3..5]), "\n";
            }
            else
            {
                print UTGO $a[0], "\t", 0, "\t", $a[2], "\t", join("\t",@a[3..5]),"\n";
            }
        }
        
        if( $a[5] eq '-')
        {
            print UTGO $a[0], "\t", $a[1], "\t", $a[2] + $upstream_of_tss_len, "\t", join("\t",@a[3..5]), "\n";
        }
    }
    close(UTGI);
    close(UTGO);
}

sub vlookem_all
{
    my $dir = shift;
    my $self = $dir and $dir = shift if ref $dir;
    my $vars_path = shift;
    my $an_vars = shift;
    
    opendir(DD,$dir) or die "no $dir: $!\n";
    my @dd = readdir(DD);
    @dd = grep{!/\.$/}@dd;
    my $cc = 1;
    
    open(COL,">collabels.txt") or die "no collabels\n";
    
    for(my $i = 0;$i < scalar(@dd);$i++)
    {
        my $cc2 = $cc + 1;
        
        print STDERR "working on $dd[$i]\n";
        print COL $dd[$i], "\n";
        
        `overlapSelect $dir/$dd[$i] $vars_path -idOutput out`; 
        $parsing->vlookup("$an_vars/tt$cc",1,"out",1,2,"y","$an_vars/tt$cc2");
        $cc++;
    }
    close(COL);
}

sub my_up_down_tss
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $upstream_of_tss_len = shift;
    my $downstream_of_tss_len = shift;
    my $outfile = shift;
    
    open(MUDTI,"$infile") or die "Can't open $infile: $!\n";
    open(MUDTO,">$outfile") or die "Can't open $outfile: $!\n";
    while(my $r = <MUDTI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        
        if( $a[5] eq '+')
        {
            if($a[1] > $upstream_of_tss_len)
            {
                #be cautious of $upstream of tss length
                print MUDTO $a[0], "\t", $a[1]-$upstream_of_tss_len, "\t", $a[2], "\t", join("\t",@a[3..5]), "\n";
            }
            else
            {
                print MUDTO $a[0], "\t", 0, "\t", $a[2],  "\t", join("\t",@a[3..5]),"\n";
            }
        }
        
        if( $a[5] eq '-')
        {
            if($a[2] > $downstream_of_tss_len)
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
    close(MUDTI);
    close(MUDTO);
}

sub mk_coding
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(MCI,"$infile") or die "Can't open $infile: $!\n";
    open(MCO,">$outfile") or die "Can't open $outfile: $!\n";
    
    while(my $r = <MCI>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        
        $a[3] =~ /^([^\-]+)\-([^\-]+)\-(.*)$/;
        
        my $tt = $3;
        
        $tt =~ s/\|/,/;
        my @w2 = split(",",$tt);
        
        if($w2[1] =~ /[+\-]/)
        {
            
            my $t = $w2[1];
            
            $t =~ s/^[\+\-]//;
            
            my $ll = length($t)/3;
            if($ll =~ /\./) # no triplet - disruptive
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
    close(MCI);
    close(MCO);
}

sub Syn
{
    my $SNP = shift;
    my $self = $SNP and $SNP = shift if ref $SNP;
    my $genome = shift;
    my $outfile = shift;
    my $dbpath = shift;
    
    open(SO,">$outfile") or die "Can't open $outfile: $!\n";
    
    my $syn = 0;
    my $nsyn = 0;
    my $ccds;
    my $ccds_tab;
    
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
    
    if($genome eq 'hg19')
    {
        $ccds = "$dbpath/ccds.hg19.bed";
        print "No $ccds is found!\n" unless(-e $ccds);
        $ccds_tab = "$dbpath/ccds.hg19.tab";
        print "No $ccds_tab is found!\n" unless(-e $ccds_tab);    
    }
    else
    {
        print STDERR "Genome is either unknown or is not used within this pipeline.\n";
        exit;
    }

    print STDERR "overlapping with ccds.bed\n";
    my $s = "overlapSelect $ccds $SNP ccds.input.out -mergeOutput";
    print $s;
    `$s`;
    
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
    
    open(F,"ccds.input.out") or die "no ccds.input.out\n";
    
    #make lookup
    print STDERR  "making sequence lookup table\n";
    my %hash = ();
    open(LOOK,$ccds_tab) or die "no ccds file\n";
    while(my $r = <LOOK>)
    {
        $r =~ /(.*?)\t(.*?)\n/;#hash of ccds id and fasta seq;
        $hash{$1} = $2;
    }
    
    print STDERR  "Going through your file\n";
    while(my $r = <F>)
    {
        my @a = split("\t",$r);
        
        # what are the two observed bases?
        my $aa;
        my $bb;
        if( $a[3] =~ /([ACGT])\|([ACGT])/)
        {
            $aa = $1;
            $bb = $2;
        }
        else
        {
            print $a[3], " does not have any base identities\n";
            next;
        }
        
        if($aa eq '0')
        {
            next;
        }
        if($a[9] eq '-')
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
        if(exists $hash{$a[7]})
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
            
            if($can_aa eq $new_aa)
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
    
    print STDERR "Total number of events: ", $nsyn + $syn, "\n";
    print STDERR  "Synonymous substitutions: $syn\n";
    print STDERR  "NonSynonymous substitutions: $nsyn\n";
    if($syn == 0 && $nsyn == 0)
    {
        print STDERR  "Synonymous/NonSynonymous: 0\n";
    }
    else
    {
       print STDERR  "Synonymous/NonSynonymous: ", sprintf "%0.*f", 2,$syn/$nsyn;
       print "\n";
    }
    
    `rm ccds.input.out`;
    close(SO);
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
    
    for(my $i = 1;$i < scalar(@ss);$i++)
    {
        if(($poss >= $ss[$i-1]) & ($poss < $ss[$i]))
        {
            #add up exon lenghts to here
            for(my $j = 0;$j < $i-1;$j++)
            {
                $txs += $ll[$j];
            }
            #add up the exon length from the start of the next exon to muation pos;	   
            $txs += ($poss - $ss[$i-1]);
            last;#stop further searching;    
        } # if in this exons
    } #end for
    
    # if orient backwards, invert start
    if($orient eq '-')
    {
        my $txlen = 0;
        for(my $i = 0;$i < scalar(@ll);$i++)
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
    if($duss eq '')
    {
        return($start,0);
    }#perfect matching the end of ABC;
    elsif($duss =~ /3/)
    {
        return($start,1);
    }#matching the nt after ABCa;
    else
    {
        return($start,2);
    }#matching the second nt after ABCab;
}

sub get_bases
{
    #Not used in Syn.pl;

    my($line) = @_;
    $line =~ s/\|/\,/g;
    my @ll = split(",",$line);
    
    my @y = ($ll[5],$ll[6],$ll[7],$ll[8]);
    @y = sort({$b <=> $a} @y);

    my $canonicalBase = $ll[0];
    my $altBase = '';
    my $HighScore_base = '';
    my $ratio = 0;
    my $dominant_base = 'Canonical';
    my $canonicalBase_score;
    
    # What is the score for the canonical base?  Use this to figure out the
    # score ratio between canonical and non-canonical

    if($canonicalBase =~ /A/) {$canonicalBase_score = $ll[5];}
    if($canonicalBase =~ /C/) {$canonicalBase_score = $ll[6];}
    if($canonicalBase =~ /G/) {$canonicalBase_score = $ll[7];}
    if($canonicalBase =~ /T/) {$canonicalBase_score = $ll[8];}
    # Print canonical base score:
    print $canonicalBase_score, "\t";

    # First determine if the highest scoring base is non-canonical base, if so,
    # return the canonical base and the highest scoring base

    if($y[0] == $ll[5]) {$HighScore_base = 'A';}
    if($y[0] == $ll[6]) {$HighScore_base = 'C';}
    if($y[0] == $ll[7]) {$HighScore_base = 'G';}
    if($y[0] == $ll[8]) {$HighScore_base = 'T';}

    unless($HighScore_base eq $canonicalBase)
    {
	$dominant_base = 'nonCanonical';
	$ratio = ($canonicalBase_score/($y[0] + $y[1]));
	return($canonicalBase,$HighScore_base,$dominant_base,$ratio,$canonicalBase_score);
    }

    # If the highest scoring base is the canonical base, find out what the next
    # highest scoring base is and return the canonical base and the next highest
    # scoring base
    if($y[1] == $ll[5])
    {
        $altBase = 'A';
    }
    if($y[1] == $ll[6])
    {
        $altBase = 'C';
    }
    if($y[1] == $ll[7])
    {
        $altBase = 'G';
    }
    if($y[1] == $ll[8])
    {
        $altBase = 'T';
    }

    
    unless($HighScore_base eq $canonicalBase)
    {
	$dominant_base = 'Canonical';
	$ratio = ($canonicalBase_score/($y[0] + $y[1]));
	return($canonicalBase,$HighScore_base);
    }
    return(0,0);
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
    
    if($can_aa eq $new_aa)
    {
        return 0;
    }
    if($new_aa eq 'Stop')
    {
        return 3;
    }
    
    my $hydrophobic = 'AILMFWYV';
    my $polar = 'STNQ';
    my $positive = 'RHK';
    my $negative = 'DE';
    
    # the rest of them are weird so a sub to these = missense 2
    
    if(($hydrophobic =~ /$can_aa/)&($hydrophobic =~ /$new_aa/))
    {
        return 1;
    }
    if(($polar =~ /$can_aa/)&($polar =~ /$new_aa/))
    {
        return 1;
    }
    if(($positive =~ /$can_aa/)&($positive =~ /$new_aa/))
    {
        return 1;
    }
    if(($negative =~ /$can_aa/)&($negative =~ /$new_aa/))
    {
        return 1;
    }
    
    return 2;
}
############################End of Syn subs#############################################

sub mk_uniq_label
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(MULI,"$infile") or die "Can't open $infile: $!\n";
    open(MULO,">$outfile") or die "Can't open $outfile: $!\n";
    while(my $r = <MULI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print MULO $a[0]."-".$a[1], "\t", $r, "\n";
    }
    close(MULI);
    close(MULO);
}

sub process_syn
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(PSI,"$infile") or die "Can't open $infile: $!\n";
    open(PSO,">$outfile") or die "Can't open $outfile: $!\n";
    while(my $r = <PSI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print PSO $a[0]."-".$a[1], "\t", $a[9].$a[12].$a[10],
        "\t", $a[11], "\n";
    }
    close(PSI);
    close(PSO);
}

sub process_shifts
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(PRSI,"$infile") or die "Can't open $infile: $!\n";
    open(PRSO,">>$outfile") or die "Can't open $outfile: $!\n";
    while(my $r = <PRSI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        $a[3] =~ s/\|/,/;
        my @w=split(",",$a[3]);
        my $t = $w[1];
        if($t eq 'frameshift')
        {
            print PRSO $a[0]."-".$a[1], "\t", 'frameshift', "\t", 3, "\n";
        }
        elsif($t eq 'triplet_shift')
        {
            print PRSO $a[0]."-".$a[1], "\t", 'triplet_shift', "\t", 2, "\n";
        }
        else
        {
            print "shleisse \n";
            next;
        }
    }
    close(PRSI);
    close(PRSO);
}

sub pt
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(PTI,"$infile") or die "Can't open $infile: $!\n";
    open(PTO,">$outfile") or die "Can't open $outfile: $!\n";
    while(my $r = <PTI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        $a[0] =~ /^([^\-]+)\-([^\-]+)\-/;
        
        my $xx=$2-1;
        
        print PTO $1."-".$xx, "\t", $r, "\n";
    }
    close(PTI);
    close(PTO);
}

1;
