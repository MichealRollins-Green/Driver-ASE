#!/usr/bin/perl

package TCGA_Lib::TCGA_ASE_Analysis;

use warnings;
use strict;
use Parsing_Routines;
use FileHandle;
use File::Copy qw(copy);
use MCE::Map;
use Cwd 'realpath';
use Cwd;
use Math::CDF qw(:all);
use FindBin qw($Bin);
use lib "$Bin/";
no strict "refs";

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

my $parsing = TCGA_Lib::Parsing_Routines->new;

###############################################################################
# TCGA_ASE_Analysis.pm version 1
# Written by Tomas Babak, Zhongshan Cheng, Michael Vermeulen and Micheal Rollins-Green
#  in the Department of Biology at Queen's University.
#
# TCGA ASE Analysis
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
    
sub compile_ase_no_tum_norm
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $cds_dir = shift;
    my $mpileups_path = shift;
    my $ase_counts = shift;
    
    $cds_dir =~ s/\/$//;
    open(my $CANTN,$infile) or die "Can't open $infile: $!.\n";
    my @TN = <$CANTN>;
    mce_map
    {
        chomp($_);
        my @a = split("\t",$_);
        my $mpileup = $a[0];
        my $cd = $a[1];
        # have both tumor and normal pileups
        print STDERR "working on $cd\n";
        pileup_at_cd("$mpileups_path/$mpileup","$cds_dir/$cd","$ase_counts/$mpileup")
    }@TN;
}

sub pileup_at_cd
{
    #get ncbi2ucsc lookup
    #my %ucsc=get_ncbi2ucsc();
    
    # load in cd
    #print STDERR "loading cd\n";
    
    #Import contents of cds_bed files
    #0: chr1            chr1
    #1: 245042304       245042311
    #2: 245042305       245042312
    #3: T|C             A|G
    my $mpileupfile = shift;
    my $cfile = shift;
    my $ase_count = shift;
    
    my %CD;
    open(MP,"$mpileupfile") or die "Can't open $mpileupfile: $!\n";
    open(CD,"$cfile") or die "no cd: $!\n";
    open(ASEC,">$ase_count") or die "Can't open $ase_count: $!\n";
    while(my $r = <CD>)
    {
        chomp($r);
        my @a = split("\t",$r);
        $CD{$a[0]."-".$a[2]} = $a[3];
    }
    close(CD);
    
    while(my $r = <MP>)
    {
        chomp($r);
        my @a = split("\t",$r);
	
        unless(exists($CD{$a[0]."-".$a[1]}))
	{
	    next;
	}
        my $cd = $CD{$a[0]."-".$a[1]};
            my($a_counts,$b_counts) = get_counts($cd,$a[4]);
            if( ($a_counts == 0) && ($b_counts == 0))
	    {
		next;
	    }    
            #Add prefix 'chr' to $a[0] if it does not have!
            #Which is important for later calculation of gene level ase.
            #As the bed file will be intersected with the annotation bed
            $a[0]="chr".$a[0] if $a[0] =~ /^\d+/i;
            $a[0]="chrX" if ( ($a[0] =~ /^x/i) or ($a[0] =~ /23/));
            print ASEC $a[0], "-", $a[1]-1, "\t", $a[0], "\t", $a[1]-1, "\t", $a[1], "\t", $cd, "\t", $a_counts, "\t", $b_counts, "\n";
    }
    close(MP);
    close(ASEC);
}

####################pileup_at_cd subs############################
sub get_counts
{
    my($alts,$string) = @_;
    
    $alts =~ s/\|/\,/;
    my @w = split(",",$alts);
    
    my $aa=0;
    my $bb=0;

    if (defined($string))
    {
        $string =~ tr/[a-z]/[A-Z]/;
        
        my @ST = split("",$string);
        
        for (my $i=0;$i < scalar(@ST);$i++)
        {
	    if($ST[$i] eq $w[0])
	    {
	        $aa++;
	    }
	    elsif($ST[$i] eq $w[1])
	    {
	        $bb++;
	    }
	}
    }
    return($aa,$bb);
}

sub get_ncbi2ucsc 
{
    my %ucsc;
    for (my $i = 1;$i < 23;$i++)
    {
        my $nn = '00000'.$i;
        $nn = substr($nn,length($nn)-6,6);    
        $ucsc{'NC_'.$nn} = $i;
    }
	
    return(%ucsc);	
}
####################end of pileup_at_cd subs############################

sub change_beds_gl_ase
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $RNA_Path_cds_sorted = shift;
    my $ase_cds_sorted_ase = shift;

    mce_map_f
    {
	my $TCGA = [split("\t",$_)]->[0];
        chomp($TCGA);
        print "Now working on $TCGA\n";
        Change_CDChr_Info_For_GeneASE("$TCGA","$RNA_Path_cds_sorted","$ase_cds_sorted_ase/$TCGA");
    }$infile;
}

sub Change_CDChr_Info_For_GeneASE
{
    my $cd = shift;
    my $cd_dir = shift;
    my $outdir = shift;
    
    $cd_dir =~ s/\/$//;
    
    my $FH = FileHandle->new("$cd_dir/$cd") or die "no $cd_dir/$cd: $!\n";
    my $OU = FileHandle->new(">$outdir") or die "Can't open $outdir: $!\n";
    while(my $l = <$FH>)
    {
        chomp($l);
        $l =~ s/^(chr23|23|X)/chrX/i;
        $l =~ s/^(\d+)/chr$1/i;
        my @as = split("\t| ",$l);
        print $OU join("\t",@as[0..3]),"\n"; #remove rsID to avoid of downstream interuption;
    }
    $FH->close;
    $OU->close;
}

sub compile_gene_ase_faster
{
    #Turn off strict refs to inable file handle contain numbers;
    #no strict "refs";
    
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $cds_ase_dir = shift;
    my $ref_bed = shift;
    my $ase_path = shift;
    my $ase_counts = shift;
    my $gene_level = shift;
    
    $cds_ase_dir =~ s/\/$//;
    
    open(my $gene,"$infile") or die "Can't open $infile: $!\n";

    my @gene_level = <$gene>;
    mce_map
    {
        chomp(my $r = $_);
        my @a = split("\t",$r);
        my $lines = join(":",@a);
        my $rr = rand();
        $rr = substr($rr,2,6); 
        print STDERR "working on $a[0]\n";
      
      #compile_gene_ase($rr,$ref_bed,$cds_dir,$r);
      #when use the above sub, be caution about pbinom from the perl Math::CDF.
      #For the following bsub, R pbinom.R was used instead of perl pbinom.
      
      compile_gene_ase("$rr","$ref_bed","$cds_ase_dir","$lines","$ase_path","$ase_counts","$gene_level");
      
    }@gene_level;
    close($gene);
}

#########################################compile_gene_faster subs#############################################
sub compile_gene_ase
{
    my $rd_id = shift;#random file id;
    my $bed_ref = shift;
    my $cds_ase_dir = shift;
    my $l = shift;
    my $ase_path = shift;
    my $ase_counts = shift;
    my $gene_level = shift;
    
    my @es = split(":",$l);

    my $temp = "temp";

    mkdir "$ase_path/$temp" unless(-d "$ase_path/$temp");
   
    #rejoins because the split above separates at the colon between the ase_counts files and cds_ase_sorted files
    #which means that it will alsp split the ase sorted files as they have a colon in them for UUID and TCGA ID.
    #$es[0] = join(".",$es[0],$es[1]);
   
    convert2mm("$ase_path/$ase_counts/$es[0]","$ase_path/$temp/mm.$rd_id.bed");
    #Overlap ase_count bed with cds_bed;
    #Be caution about the bed format between the two files, especially for chr label!
    TallyOverDefBed("$ase_path/$temp/mm.$rd_id.bed","$cds_ase_dir/$es[1]","$ase_path/$temp/tallied.$rd_id.bed","$ase_path/$temp");
    AlleleCount2bed("$ase_path/$temp/tallied.$rd_id.bed","$bed_ref","$ase_path/$gene_level/$es[0]","$ase_path/$temp");   
}

sub convert2mm
{
    #summarize reads of all four bases (A,C,G,T)and prepare data for TallyOVerDefBed
    my $ase_count_file = shift;
    my $out_file = shift;
    
    my $FH = "FH".rand();
    my $OUT = "OUT".rand();
    open($FH, "$ase_count_file") or die "Can not open the ase count file $ase_count_file $!\n";
    open($OUT, ">$out_file") or die "Can not write data into the file $out_file  $!\n";
    
    my @Newmm = map
    {
        my $r = $_;
        chomp($r);
        my @a = split("\t",$r);
     
        my $A = 0;
        my $C = 0;
        my $G = 0;
        my $T = 0;
     
        $a[4] =~ s/\|/,/;
        my @w = split(",",$a[4]);
     
        my $aa = $w[0];
        my $bb = $w[1];
     
        eval('$'.$aa.'+='.$a[5]);
        eval('$'.$bb.'+='.$a[6]);
     
        my $string = $A."\|".$C."\|".$G."\|".$T;
      
        my $out = join("\t",@a[1..3])."\t"."N|".$string."\|".$string."\|".$string;
      
        print $OUT $out,"\n";
    }(<$FH>);
   close $FH;
   close $OUT;
}

sub TallyOverDefBed
{
    my $snp = shift;#snp.bed
    my $def = shift;#def.bed, such as refseq.ucsc.ensembl.mrna.hg9.nr.bed
    my $out = shift;#output filename
    my $ase_path = shift;#ase path
    
    my $sr = 0;
    my $minC = 0;
    my $phred = 0;
    my $minUC = 0;
    my $minUCnr = 0;
   
    # random filename to allow multi-jobs in same directory
    my $rr = rand();
    my $rr2 = rand();
    $rr = substr($rr,2,6).substr($rr2,2,6);
   
    my $OUT_FH = "OUT_FH".rand();
    open($OUT_FH,">$out") or die "Can not open the file to write data: $!";

    unless ($phred =~ /^[0-9\.]+$/)
    {
        print "-phred must be a number, not $phred\n";
        exit;
    }
    unless ($minC =~ /^[0-9]+$/)
    {
        print "-minC must be a number, not $phred\n";
        exit;
    }
    
    if($sr == 1)
    {
        #split file
        my $PFH = "P".rand();
        my $MFH = "M".rand();
        my $FFH = "F".rand();
       
        open($PFH,">$ase_path/plus.$rr.bed") or die "cant write data into plus.$rr.bed: $!";;
        open($MFH,">$ase_path/minus.$rr.bed") or die "cant write data into minus.$rr.bed: $!";
        open($FFH,$snp) or die "$snp not here: $!\n";
        while(my $r = <$FFH>)
        {
            if($r =~ /\t\+[\s]/)
            {
                print $PFH $r;
            }
            else
            {
                print $MFH $r;
            }
        }
       	close $FFH;
       	close $MFH;
        close $PFH;
    }
    
    # overlap with def.bed
    
    if($sr == 1)
    {
        `overlapSelect $def $ase_path/plus.$rr.bed -mergeOutput $ase_path/tmp.$rr.plus.out`;
        `overlapSelect $def $ase_path/minus.$rr.bed -mergeOutput $ase_path/tmp.$rr.minus.out`;
        my $F_FH = "F_FH".rand();
        open($F_FH,"$ase_path/tmp.$rr.plus.out") or die "no tmp.plus.out: $!\n";
        
        while(my $r = <$F_FH>)
        {    
            chomp($r);
            my @aa = split("\t",$r);
            my $x = $aa[3];
            $x =~ s/\|/,/g;
            my @bb = split(",",$x);
            my $A = $bb[1];
            my $C = $bb[2];
            my $G = $bb[3];
            my $T = $bb[4];
            my $qA = $bb[5];
            my $qC = $bb[6];
            my $qG = $bb[7];
            my $qT = $bb[8];
            my $cA = $bb[9];
            my $cC = $bb[10];
            my $cG = $bb[11];
            my $cT = $bb[12];
       	    
            my $y = $aa[9];
            $y =~ s/\|/,/g;
	    
	    if ($y =~ /-/)
	    {
		next;
	    }
            
	    my @perl = split(",",$y);
           
            my $can = eval("\$$perl[0]");
            my $alt = eval("\$$perl[1]");
       	    my $qcan = eval("\$q$perl[0]");
            my $qalt = eval("\$q$perl[1]");
            my $ccan = eval("\$c$perl[0]");
            my $calt = eval("\$c$perl[1]");
           
       	    if($phred > 0)
            {
                #want a phred cut
       		
           	my($nonDB1,$nonDB2) = getNonDB($perl[0],$perl[1]);
               	my $nonDB1q = eval("\$q$nonDB1");
           	my $nonDB2q = eval("\$q$nonDB2");
           	unless(($nonDB1q >= $phred)|($nonDB2q >= $phred))
                {
       		    unless(($can >= $minC) & ($alt >= $minC))
                    {
                        next;
                    }								 
                    unless(($calt >= $minUC) & ($ccan >= $minUC))
                    {
                        next;
                    }
       		    unless($calt >= $minUCnr)
                    {
                        next;
                    }
                    #print chrnum,st,end,transcriptID, and read counts;			   
                    print $OUT_FH join("\t",@aa[0..2]), "\t", $aa[9]."\|".$can."\|".$alt."\|".$qcan."\|".$qalt."\|".$ccan."\|".$calt, "\t0\t+\n";
               	}
       			
       	    } #if
       	    else
            {
                # no phred cuts
       		unless(($can >= $minC) & ($alt >= $minC))
                {
                    next;
                }	
       		unless(($calt >= $minUC) & ($ccan >= $minUC))
                {
                    next;
                }
       		unless($calt >= $minUCnr)
                {
                    next;
                }
		#print chrnum,st,end,transcriptID, and read counts;			
               	print $OUT_FH join("\t",@aa[0..2]), "\t", $aa[9]."\|".$can."\|".$alt."\|".$qcan."\|".$qalt."\|".$ccan."\|".$calt, "\t0\t+\n";
       	    }
        } #while plus
        close $F_FH;
        `rm $ase_path/tmp.$rr.plus.out`;
           
        my $FH = "FH".rand();
        open($FH,"$ase_path/tmp.$rr.minus.out") or die "no tmp.minus.out: $!\n";
           
        while(my $r = <$FH>)
        {    
            chomp($r);
            my @aa = split("\t",$r);
            my $x = $aa[3];
            $x =~ s/\|/,/g;
            my @bb = split(",",$x);
            my $A = $bb[1];
            my $C = $bb[2];
            my $G = $bb[3];
            my $T = $bb[4];
            my $qA = $bb[5];
            my $qC = $bb[6];
            my $qG = $bb[7];
            my $qT = $bb[8];
            my $cA = $bb[9];
            my $cC = $bb[10];
            my $cG = $bb[11];
            my $cT = $bb[12];
            
            my $y = $aa[9];
            $y =~ s/\|/,/g;
	    
	    if ($y =~ /-/)
	    {
		next;
	    }
            
	    my @perl = split(",",$y);
           
            my $can = eval("\$$perl[0]");
            my $alt = eval("\$$perl[1]");
            my $qcan = eval("\$q$perl[0]");
            my $qalt = eval("\$q$perl[1]");
            my $ccan = eval("\$c$perl[0]");
            my $calt = eval("\$c$perl[1]");
               	
            if($phred > 0)
            {
                #want a phred cut
       		
           	my($nonDB1,$nonDB2) = getNonDB($perl[0],$perl[1]);
               	my $nonDB1q = eval("\$q$nonDB1");
           	my $nonDB2q = eval("\$q$nonDB2");
           	unless(($nonDB1q >= $phred)|($nonDB2q >= $phred))
                {
                    unless(($can >= $minC) & ($alt >= $minC))
                    {
                        next;
                    }	
       		    unless(($calt >= $minUC) & ($ccan >= $minUC))
                    {
                        next;
                    }
       		    unless($calt >= $minUCnr)
                    {
                        next;
                    }
                    #print chrnum,st,end,transcriptID, and read counts;			   
                    print $OUT_FH join("\t",@aa[0..2]), "\t", $aa[9]."\|".$can."\|".$alt."\|".$qcan."\|".$qalt."\|".$ccan."\|".$calt, "\t0\t-\n";
               	}
       	    } #if
       	    else
            {
                # no phred cuts
       		unless(($can >= $minC) & ($alt >= $minC))
                {
                    next;
                }	
       		unless(($calt >= $minUC) & ($ccan >= $minUC))
                {
                    next;
                }
       		unless($calt >= $minUCnr)
                {
                    next;
                }
		#print chrnum,st,end,transcriptID, and read counts;			
               	print $OUT_FH join("\t",@aa[0..2]), "\t", $aa[9]."\|".$can."\|".$alt."\|".$qcan."\|".$qalt."\|".$ccan."\|".$calt, "\t0\t-\n";
       	    }
        }
       	close $FH;
       	unlink("$ase_path/tmp.$rr.minus.out");
    }
    else
    {
        # no -strand
       	`overlapSelect $def $snp -mergeOutput $ase_path/tmp.$rr.all.out`;
        my $FH_= "FH_".rand();
        open($FH_,"$ase_path/tmp.$rr.all.out") or die "no tmp.$rr.all.out: $!\n";
           
        while(my $r = <$FH_>)
        {    
            chomp($r);
            my @aa = split("\t",$r);
            my $x = $aa[3];
            $x =~ s/\|/,/g;
            my @bb = split(",",$x);
            my $A = $bb[1];
            my $C = $bb[2];
            my $G = $bb[3];
            my $T = $bb[4];
            my $qA = $bb[5];
            my $qC = $bb[6];
            my $qG = $bb[7];
            my $qT = $bb[8];
            my $cA = $bb[9];
            my $cC = $bb[10];
            my $cG = $bb[11];
            my $cT = $bb[12];
            
            my $y = $aa[scalar(@aa)-1];
            $y =~ s/\|/,/g;
	    
	    if ($y =~ /-/)
	    {
		next;
	    }
	    
            my @perl = split(",",$y);
           
            my $can = eval("\$$perl[0]");
            my $alt = eval("\$$perl[1]");
            my $qcan = eval("\$q".$perl[0]);
            my $qalt = eval("\$q".$perl[1]);
            my $ccan = eval("\$c$perl[0]");
            my $calt = eval("\$c$perl[1]");
	    
       	    if($phred > 0)
            {
                #want a phred cut
       		
           	my($nonDB1,$nonDB2) = getNonDB($perl[0],$perl[1]);
               	my $nonDB1q = eval("\$q$nonDB1");
           	my $nonDB2q = eval("\$q $nonDB2");
           	unless(($nonDB1q>=$phred)|($nonDB2q>=$phred))
                {
               	    unless(($can >= $minC) & ($alt >= $minC))
                    {
                        next;
                    }	
       		    unless(($calt >= $minUC) & ($ccan >= $minUC))
                    {
                        next;
                    }
       		    unless($calt >= $minUCnr)
                    {
                        next;
                    }
		    #print chrnum,st,end,transcriptID, and read counts;
               	    print $OUT_FH join("\t",@aa[0..2]), "\t", $aa[7]."\|".$can."\|".$alt."\|".$qcan."\|".$qalt."\|".$ccan."\|".$calt."\n";
               	}
       		
       	    } #if
       	    else
            {
                # no phred cuts
       		unless(($can >= $minC) & ($alt >= $minC))
                {
                    next;
                }	
       		unless(($calt >= $minUC) & ($ccan >= $minUC))
                {
                    next;
                }
       		unless($calt >= $minUCnr)
                {
                    next;
                }
		#print chrnum,st,end,transcriptID, and read counts;			
               	print $OUT_FH join("\t",@aa[0..2]), "\t", $aa[7]."\|".$can."\|".$alt."\|".$qcan."\|".$qalt."\|".$ccan."\|".$calt."\n";
       	    }
       	    
        }
        close $FH_;
       	unlink("$ase_path/tmp.$rr.all.out");
    }
}

sub AlleleCount2bed
{   
    my $talliedBed = shift;
    my $RefBed = shift;
    my $gene_ase_out = shift;
    my $ase_path = shift;
    
    my $newRefBed = $RefBed;
    my $rn = rand();
    $rn = substr($rn,2,6);
    $newRefBed =~ s/\.bed/$rn\.bed/;
    system("cp $RefBed $newRefBed");
    $RefBed = $newRefBed;
   
    my $rr = rand();
    my $rr2 = rand();
    $rr = substr($rr,2,6).substr($rr2,2,6);
    #########################################################################################
    # computes counts over transcripts from                                                 #
    #tallied.bed combineAcrossCoord.bed {options}                                           #
    #                                                                                       #
    #  options:                                                                             #
    #                                                                                       #
    #  -selectRange    will consider entire range in combineAcrossCoord.bed entries         #
    #  -out filename   default=counts.bed                                                   #
    #  -strand         only works if both alleleCounts.be and coord.bed have strand info    #
    #  -snps           prints the raw SNP counts to one column/gene                         #
    #  -agree          prints proportion of SNPs with 5+ counts in agreement on direction   #
    #                  and number of informative snps (5+ counts)                           #
    #  -highest_b       reports snp and p-value of highest-b snp                            #
    #  -snp_ids        also prints snp ids in same sequence as snps                         #
    #  -print_cd_snp   print A|C,G|T,C|C -> het ids used                                    #
    #                                                                                       #
    #  N.B. combineAcrossCoord.bed must have Unique IDs in column 4                         #
    #########################################################################################
    
    my $old_print = '';       
    my $strand = '';#'-strand'
    my $sr = '-selectRange';#'-selectRange'
    my $out = $gene_ase_out;
    my $snps = 'y';#'y'
    my $agree = 'y';#'y'
    my $highest_b = 'n';#'y'
    my $snp_ids = 'y';#'y'
    my $print_cd = 'y';#'y'
    
    `overlapSelect $RefBed $talliedBed $sr $strand -mergeOutput $ase_path/temp.$rr.bed`;

    #Contents of temp.$rr.bed;
    #Array_Num: Headers         RowEelments
    #0: chr6                    chr6
    #1: 116590801               116590806
    #2: 116590802               116590807
    #3: T|A|4|0|4|0|4|0         T|G|4|0|4|0|4|0
    #4: chr6                    chr6
    #5: 116590506               116590506
    #6: 116592627               116592627
    #7: AB075478                AB075478
    #8: 0                       0
    #9: -                       -
    #10: 116590506              116590506
    #11: 116592627              116592627
    #12: 0                      0
    #13: 4                      4
    #14: 337,306,210,1266,      337,306,210,1266,
    #15: 0,338,644,855,         0,338,644,855,
    
    my $g = `head -n 1 $talliedBed`;
    my @g = split("\t",$g);
    my $ff = scalar(@g);
    #conteins of cloumns 1-4: chr1    200295519       200296221        T|C|100|182|100|182|100|182
    $ff += 4;
    #The following SORT step is necessary for correct calculation of parameters, such as ppp,logr, and so on;
    `sort -k $ff,$ff $ase_path/temp.$rr.bed > $ase_path/sorted.$rr.temp.bed`;
    `mv $ase_path/sorted.$rr.temp.bed $ase_path/temp.$rr.bed`;
    
    my $BEDFH = "BEDFH".rand();
    open($BEDFH,"$ase_path/temp.$rr.bed") or die "Cannot open $ase_path/temp.$rr.bed: $!\n";
    my $OUTFH = "OUTFH".rand();
    open($OUTFH,">$out") or die "Cannot open $out for writing: $!\n";
    
    my $r = <$BEDFH>;
    $r =~ s/\n//g;
    
    my @x = split("\t",$r);
    my $old = $x[$ff-1];#transcript Id;
    my $chr = $x[0];
    my $start = $x[$ff-3];
    my $end = $x[$ff-2];
    my $orient = $x[$ff+1];
    my @pr = @x;
    
    $x[3] =~ s/\|/,/g;#T|C|100|182|100|182|100|182
    my @tt1;
    my $tt2;
    my @counts = split(",",$x[3]);
    my $can = $counts[2];#100 of T|C|100|182|100|182|100|182
    my $var = $counts[3];#182 of T|C|100|182|100|182|100|182
    my $snps_out = $counts[2]."|".$counts[3];#100|182 of T|C|100|182|100|182|100|182
    my $snp_ids_out = $x[0]."-".$x[1];#chrid-snp_st;
    $tt2 = $x[3]; $tt2=~s/\|/,/g; @tt1=split(",",$tt2); my $print_cd_out=$tt1[0]."|".$tt1[1];#allele1|allele2
    my $aa_dir = 0;
    my $bb_dir = 0;
    my $informative_snps = 0;
    
    unless($counts[2] == $counts[3])
    {
        
	if(($counts[2] + $counts[3]) > 4)
        {
	    $informative_snps++;
	    if($counts[2] > $counts[3])
            {
	        $aa_dir++;
	    }
	    else
            {
		$bb_dir++;
	    }
	}
    }
    
    while($r = <$BEDFH>)
    {
        chomp($r);
        @x = split("\t",$r);
        
        my $incr = $x[$ff-1];#transcript id;
        
        if ($incr eq $old)
        {
            #compare with old transcript id;           
            $x[3] =~ s/\|/,/g;
            @counts = split(",",$x[3]);
            $can += $counts[2];
            $var += $counts[3];
            $snps_out = $snps_out.",".$counts[2]."|".$counts[3];
            $snp_ids_out = $snp_ids_out.",".$x[0]."-".$x[1];
            $tt2 = $x[3];
	    $tt2=~s/\|/,/g;
	    @tt1=split(",",$tt2);
	    $print_cd_out = $print_cd_out.",".$tt1[0]."|".$tt1[1];
            unless($counts[2] == $counts[3])
            {
		if(($counts[2] + $counts[3]) > 4)
                {
		    $informative_snps++;
		    if($counts[2] > $counts[3])
                    {
			$aa_dir++;
		    }
		    else
                    {
			$bb_dir++;
		    }
		}
	    }
        }   
        else
        {
            my $ppp = get_p($can,$var);
            if($can < $var)
            {
		$ppp = $ppp*(-1);
	    }	
            my $logr = 0;
            unless(($can == 0)|($var == 0))
            {
		$logr = log($can/$var)/log(2);
	    }
            print $OUTFH $chr, "\t", $start, "\t", $end, "\t", $old, "\t", $can,
	    "\t", $var, "\t", $ppp, "\t", $logr, "\t", $orient;
            
            #if corodBed file has exon coordinates, print them
            if($snps eq 'y')
            {
                print $OUTFH "\t", $snps_out;
            }
            
            if($agree eq 'y')
            {
                my $agree_out = agree($aa_dir,$bb_dir);
                my $num_snps=num($snps_out);
                print $OUTFH "\t", $agree_out, "\t", $informative_snps;
	    }
            
            if($highest_b eq 'y')
            {
		print $OUTFH "\t", highest_b($snps_out,"$ase_path");
	    }
            
            if($snp_ids eq 'y')
            {
		print $OUTFH "\t", $snp_ids_out;
	    }
            
            if($print_cd eq 'y')
            {
		print $OUTFH "\t", $print_cd_out;
	    }
            
            if( (scalar(@pr) <= 18)|($sr eq '-selectRange'))
            {
		print $OUTFH "\n";
	    }
            else
            {
		print $OUTFH "\t", join("\t",@pr[18..(scalar(@pr)-1)]), "\n";
	    }
            
	    #To process data of new transcript id;
            $x[3] =~ s/\|/,/g;
            @counts = split(",",$x[3]);
            $can = $counts[2];
            $var = $counts[3];
            $snps_out = $counts[2]."|".$counts[3];
            $snp_ids_out = $x[0]."-".$x[1];
            $tt2 = $x[3];
	    $tt2=~s/\|/,/g;
	    @tt1=split(",",$tt2);
	    $print_cd_out=$tt1[0]."|".$tt1[1];
            $aa_dir = 0;
            $bb_dir = 0;
            $informative_snps = 0;
    	    
            unless($counts[2] == $counts[3])
            {
		if(($counts[2] + $counts[3]) > 4)
                {
		    $informative_snps++;
		    if($counts[2] > $counts[3])
                    {
			$aa_dir++;
		    }
		    else
                    {
			$bb_dir++;
		    }
		}
	    }       		
        }#end of else;
       
	#Assign new value for each variable
	#N.B.: This is very important!
        $chr = $x[0];
        $old = $incr;
        $start = $x[$ff-3];
        $end = $x[$ff-2];
        $orient = $x[$ff+1];
        @pr = @x;
        
    }#end of each line of while loop;
   
    # print last entry
    my $ppp = get_p($can,$var);
    if($can < $var)
    {
        $ppp = $ppp*(-1);
    }	
   
    my $logr = 0;
    unless(($can == 0)|($var == 0))
    {
        $logr = log($can/$var)/log(2);
    }
    
    print $OUTFH $chr, "\t", $start, "\t", $end, "\t", $old, "\t", $can, "\t", $var, "\t", $ppp, "\t", $logr, "\t", $orient;

    #if corodBed file has exon coordinates, print them
    if($snps eq 'y')
    {
        print $OUTFH "\t", $snps_out;
    }
    
    if($agree eq 'y')
    {
	my $agree_out = agree($aa_dir,$bb_dir);
	my $num_snps = num($snps_out);
	print $OUTFH "\t", $agree_out, "\t", $informative_snps;
    }

    if($highest_b eq 'y')
    {
	print $OUTFH "\t", highest_b($snps_out,"$ase_path");
    }

    if($snp_ids eq 'y')
    {
	print $OUTFH "\t", $snp_ids_out;
    }
    
    if($print_cd eq 'y')
    {
        print $OUTFH "\t", $print_cd_out;
    }

    if( (scalar(@pr) <= 18)|($sr eq '-selectRange'))
    {
	print $OUTFH "\n";
    }
    else
    {
        print $OUTFH "\t", join("\t",@pr[18..(scalar(@pr)-1)]), "\n";
    }
   
    close($OUTFH);
    close($BEDFH);
    `rm $ase_path/temp.$rr.bed`;
    system("rm $RefBed");
}

##############
# Small SUBS
##############

sub getNonDB
{  
    my($xx,$yy) = @_;

    my $cc = 'ACGT';
    $cc =~ s/$xx//i;
    $cc =~ s/$yy//i;
    return(substr($cc,0,1),substr($cc,1,1));
}


sub getC
{
    my($a) = @_;

    unless ($a =~ /^[0-9]+$/)
    {
        return 1;
    }
    else
    {
        return $a;
    }
}

sub get_p
{
    #Try to detect Perl module Math::CDF;
    BEGIN
    {
	eval 'use Math::CDF qw(:all); 1';
    }
    
    my($s,$t,$RNA_ASE) = @_;
    
    my $total = $s + $t;
    my $out = 0;
    
    if( ($s == 0) & ($t == 0))
    {
        return 0;
    }
    elsif($s < $t)
    {
	if (defined &pbinom)
        {
	    $out = pbinom($s,($s+$t),0.5);
	}
        else
        {
	    $out = `R --no-save $s $total 0.5 <$RNA_ASE/pbinom.R --quiet|tail -2|head   -1`;
	    $out =~ s/\[1\]\s*//;		
	}
    }
    else
    {
	if (defined &pbinom)
        {
	    $out = pbinom($t,($s+$t),0.5);
	}
        else
        {
	    $out = `R --no-save $t $total 0.5 <$RNA_ASE/pbinom.R --quiet|tail -2|head -1`;
	    $out =~ s/\[1\]\s*//;		
	}  
        
    }
    
    if($out == 0)
    {
        return 500;
    }
    else
    {
        return(-1*log($out)/log(10));
    }		
}

sub agree 
{
    my($aa,$bb) = @_;
    
    if(($aa+$bb) == 0)
    {
        return 0;
    }
    return(substr(max($aa,$bb)/($aa+$bb),0,6));
}

sub max
{
    my($aa,$bb) = @_;

    if($aa > $bb)
    {
        return($aa);
    }
    else
    {
        return($bb);
    }
}

sub min
{
    my($aa,$bb) = @_;

    if($aa > $bb)
    {
        return($bb);
    }
    else
    {
        return($aa);
    }
}

sub num
{
    my($in) = @_;

    my @w = split(",",$in);
    return(scalar(@w));
}

sub highest_b 
{
    my($in,$RNA_ASE) = @_;
    
    chomp($in);
    my @w = split(",",$in);
    my $out = $w[0];
    $out =~ /^([0-9]+)\|([0-9]+)/;
    
    my $b = min($1,$2);
    my $p_out = get_p($1,$2);
    if($1 < $2)
    {
        $p_out = $p_out*-1;
    }
    
    $out =~ s/\|/\t/;
    for(my $i = 1;$i < scalar(@w);$i++)
    {
        $w[$i] =~ /^([0-9]+)\|([0-9]+)/;
        my $ss = $1;
        my $tt = $2;
        $out =~ s/\|/\t/;
	my $mm = min($ss,$tt);
        if($mm > $b)
        {
            $out = $ss."\t".$tt;
	    $b = $mm;
	    $p_out=get_p($ss,$tt,$RNA_ASE);
	    if($ss < $tt)
            {
		$p_out = $p_out*-1;
	    }
	}
    }    
    return($out."\t".$p_out);
}
#########################################end of compile_gene_faster subs#############################################

sub pull_TN
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $bad_cnv_dir = shift;
    my $bad_snp_bed_dir = shift;
    my $outfile = shift;
    
    open(PTNI,"$infile") or die "Can't open $infile: $!\n";
    open(PTNO,">$outfile") or die "Can't open $outfile: $!\n";
    
    my $r = <PTNI>;
    chomp($r);
    my @old = split("\t",$r);
    my @a;
    my $bad_snp = 0;
    my $bad_cnv = 0;
    
    if(-e $bad_cnv_dir."/".$old[1])
    {
        $bad_cnv = 1;
    }
    if(-e $bad_snp_bed_dir."/".$old[1])
    {
        $bad_snp = 1;
    }
    
    while($r = <PTNI>)
    {
        chomp($r);
        @a = split("\t",$r);
        
        if ($old[3] < 10 || $old[3] < 11)
        {
            if ($old[3]!~/0/)
            {
                $old[3] = '0' . $old[3];
            }
        }
        
        if (($old[2] eq $a[2]) & ($old[3] eq '01') & ( ($a[3] eq '10') | ($a[3] eq '11'))) # have matched N
        {
            unless(($old[3] eq '10') | ($old[3] eq '11') | ($old[3] eq '06'))
            {
                print PTNO $old[0], "\t", $a[0], "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
            }
        }
        else
        {
            unless(($old[3] eq '10') | ($old[3] eq '11') | ($old[3] eq '06'))
            {
                print PTNO $old[0], "\t", "NaN", "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
            }
        }
        undef @old;
        @old = @a;
        
       if(-e $bad_cnv_dir."/".$old[1])
        {
            $bad_cnv = 1;
        }
        else
        {
            $bad_cnv = 0;
        }
        if(-e $bad_snp_bed_dir."/".$old[1])
        {
            $bad_snp = 1;
        }
        else
        {
            $bad_snp = 0;
        }
    }
    
    if(($old[3] eq '10') | ($old[3] eq '11') | ($old[3] eq '06'))
    {
	EXIT_IF:
	{
	    last EXIT_IF;
	}
    }
    if(($old[2] eq $a[2]) & ($old[3] eq '01') & (($a[3] eq '10') | ($a[3] eq '11'))) # have matched N
    {
        print PTNO $old[0], "\t", $a[0], "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
    }
    else
    {
        print PTNO $old[0], "\t", "NaN", "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
    }
    close(PTNI);
    close(PTNO);
}

sub have_cnv_snp
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(HCSI,"$infile") or die "Can't open $infile: $!\n";
    open(HCSO,">$outfile") or die "Can't open $outfile: $!\n";
    
    while(my $r = <HCSI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        if( ($a[3] == 1) && ($a[4] == 1))
        {
            print HCSO $r, "\n";
        }
    }
    close(HCSI);
    close(HCSO);
}

sub mk_files
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $ase_path = shift;
    my $gene_level = shift;
    
    open(MKFI,"$infile") or die "Can't open $infile: $!\n";
    open(N,">normal_ff") or die "no normal out\n";
    open(T,">tumor_ff") or die "no tumor out\n";
    
    while(my $r = <MKFI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        
        print T "$ase_path/$gene_level/".$a[0], "\t", substr($a[2],0,12), "\n";
        unless($a[1] eq 'NaN')
        {
            print N "$ase_path/$gene_level/".$a[1], "\t", substr($a[2],0,12), "\n";
        }
    }
    close(MKFI);
    close(N);
    close(T);
}

sub run_problem_cnvs
{
    my $bad_cnv_dir = shift;
    my $self = $bad_cnv_dir and $bad_cnv_dir = shift if ref $bad_cnv_dir;
    my $infile = shift;
    my $ase_path = shift;
    my $matrix = shift;
    my $problem_cnvs = shift;
    
    open(RPCI,"$infile") or die "Can't open $infile: $!";
    
    mkdir "$ase_path/$problem_cnvs" unless(-d "$ase_path/$problem_cnvs");
    system("rm -f $ase_path/$problem_cnvs/*");
    while(my $r = <RPCI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        
        print STDERR "working on $a[0]\n";
        
        my $bed = "$a[2].bed";
        
        `overlapSelect $bad_cnv_dir/$bed $ase_path/$matrix/rowlabels.bed $ase_path/cnv_t`;
        pull_problems("$ase_path/cnv_t","$ase_path/$problem_cnvs/$a[0]");
    }
    close(RPCI);
}

sub run_problem_snps
{
    my $bad_snp_dir = shift;
    my $self = $bad_snp_dir and $bad_snp_dir = shift if ref $bad_snp_dir;
    my $infile = shift;
    my $ase_path = shift;
    my $matrix = shift;
    my $problem_snps = shift;
    
    $bad_snp_dir =~ s/\/$//;
    
    mkdir "$ase_path/$problem_snps" unless(-d "$ase_path/$problem_snps");
    system("rm -f $ase_path/$problem_snps/*");
    
    open(PPSI,"$infile") or die "Can't open $infile: $!\n";
    
    while(my $r = <PPSI>)
    {
	chomp($r);
	my @a = split("\t",$r);
	
	print STDERR "working on $a[0]\n";
	
	my $bed = "$a[2].bed";
	
	`overlapSelect $bad_snp_dir/$bed $ase_path/$matrix/rowlabels.bed $ase_path/snp_t`;
	pull_problems("$ase_path/snp_t","$ase_path/$problem_snps/$a[0]");
    }
    close(PPSI);
}

sub pull_problems
{
    my $infile = shift;
    my $outdir = shift;
    
    open(PPI,"$infile") or die "Can't open $infile: $!\n";
    open(PPO,">$outdir") or die "Can't open $outdir: $!\n";
    while(my $r = <PPI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print PPO $a[3], "\t", 1, "\n";
    }
    close(PPI);
    close(PPO);
}

1;
