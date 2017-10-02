package Driver_ASE_Lib::Gene_Level_ASE_Analysis;

use warnings;
use strict;
use FileHandle;
use File::Copy qw(copy);
use MCE::Map;
use Cwd 'realpath';
use Cwd;
use Math::CDF qw(:all);
use FindBin qw($Bin);
use lib "$Bin/";
use Parsing_Routines;
no strict "refs";
use autodie;
no warnings 'once';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

my $parsing = Driver_ASE_Lib::Parsing_Routines->new;

###############################################################################
# Gene_Level_ASE_Analysis.pm version 1
# Written by Tomas Babak, Zhongshan Cheng, Michael Vermeulen and Micheal Rollins-Green
#  in the Department of Biology at Queen's University.
#
# Gene_Level_ASE_Analysis
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
    my $mpileup_file = shift; #file with list of mpileups files
    my $self = $mpileup_file and $mpileup_file = shift if ref $mpileup_file;
    my ($cds_dir,$mpileups_path,$ase_counts) = @_; #path to directory with sorted bed files, path to directory that hold mpileups, path to directory that will hold ase counts
    
    $cds_dir =~ s/\/$//;

    open (my $CANTN,$mpileup_file);
    my @TN = <$CANTN>;
    mce_map
    {
        chomp($_);
        my @a = split("\t",$_);
        my $mpileup = $a[0];
        my $cd = $a[2];
	
	if (-s "$mpileups_path/$mpileup" != 0)
	{
	    print "Details for ase analysis: $_\n";
	    # have both tumor and normal pileups
	    print "Working on $cd\n";
	    pileup_at_cd("$mpileups_path/$mpileup","$cds_dir/$cd","$ase_counts/$mpileup")
	}
	else
	{
	    print "$mpileups_path/$mpileup is empty!\n";
	}
    }@TN;
    close ($CANTN);
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
    my ($mpileupfile, $cd_sorted_file,$ase_count) = @_; #mpileup file, sorted bed file, directory that will hold ase count files
    
    my %CD;
    open (MP,"$mpileupfile");
    open (CD,"$cd_sorted_file");
    open (ASEC,">$ase_count");
    while(my $r = <CD>)
    {
        chomp($r);
        my @a = split("\t",$r);
        $CD{$a[0]."-".$a[2]} = $a[3];
    }
    close (CD);
    
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
            #Add prefix 'chr' to $a[0] if it does not have it!
            #Which is important for later calculation of gene level ase.
            #As the bed file will be intersected with the annotation bed
            $a[0]="chr".$a[0] if $a[0] =~ /^\d+/i;
            $a[0]="chrX" if ( ($a[0] =~ /^x/i) or ($a[0] =~ /23/));
            print ASEC $a[0], "-", $a[1]-1, "\t", $a[0], "\t", $a[1]-1, "\t", $a[1], "\t", $cd, "\t", $a_counts, "\t", $b_counts, "\n";
    }
    close (MP);
    close (ASEC);
}

####################pileup_at_cd subs############################
sub get_counts
{
    my($alts,$string) = @_;
    
    $alts =~ s/\|/\,/;
    my @w = split(",",$alts);
    
    my $aa = 0;
    my $bb = 0;

    if (defined($string))
    {
        $string =~ tr/[a-z]/[A-Z]/;
        
        my @ST = split("",$string);
        
        for (my $i = 0;$i < scalar(@ST);$i++)
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

sub compile_gene_ase_faster
{
    #Turn off strict refs to inable file handle contain numbers;
    #no strict "refs";
    
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my ($cds_dir,$ref_bed,$ase_path,$ase_counts,$gene_level,$overlapSelect) = @_; #path to sorted beds directory, path to reference bed file, path to the directory that will hold the ase analysis data, directory that holds ase count files, directory that will hold gene level analysis data, overlapSelect command
    
    $cds_dir =~ s/\/$//;
    
    open (my $gene,"$infile");
    my @gene_level = <$gene>;
    close ($gene);
    
    my $temp = "temp";
    mkdir "$ase_path/$temp" unless(-d "$ase_path/$temp");
    mce_map
    {
        chomp(my $r = $_);
        my @a = split("\t",$r);
	if (-e "$ase_path/$ase_counts/$a[0]")
	{
	    
	    print "Working on $a[0]\n";
	    my $rd_id = rand();
	    $rd_id = substr($rd_id,2,6);
	    
	    my $bed_ref = $ref_bed;
	    
	    print "Going to convert the file into mm format: $ase_path/ase_counts/$a[0]\n";
	    convert2mm("$ase_path/$ase_counts/$a[0]","$ase_path/$temp/mm.$rd_id.bed");
	    print "Failed to convert the file into mm format!\n" and next unless (-s "$ase_path/$temp/mm.$rd_id.bed" > 0);
	    
	    #Overlap ase_count bed with cds_bed;
	    #Be aware of the bed format between the two files, especially for chr label!
	    #Get rid of last column;
	    
	    Make4colBed("$cds_dir/$a[2]");
	    TallyOverDefBed("$ase_path/$temp/mm.$rd_id.bed","$cds_dir/$a[2]","$ase_path/$temp/tallied.$rd_id.bed","$ase_path/$temp","$overlapSelect");
	    AlleleCount2bed("$ase_path/$temp/tallied.$rd_id.bed","$bed_ref","$ase_path/$gene_level/$a[0]","$ase_path/$temp","$ase_path","$overlapSelect");  
	}
	else
	{
	    print "No file named $ase_path/$ase_counts/$a[0]!\n This file must be empty in mpileups directory.\n";
	}
    }@gene_level;
}

sub pull_TN
{
    my $ase_file = shift; #file with ase files listed
    my $self = $ase_file and $ase_file = shift if ref $ase_file;
    my ($bad_cnv_dir,$bad_snp_bed_dir,$cnv_snps_file,$cancer_type) = @_; #directory where bad cnvs files are stored, directory where bad snp bed files are stored, file to output results, cancer type
    
    open (PTNI,"$ase_file");
    open (PTNO,">$cnv_snps_file");
    
    my $r = <PTNI>;
    chomp($r);
    my @old = split("\t",$r);
    my @a;
    my $bad_snp = 0;
    my $bad_cnv = 0;
    
    if(-e $bad_cnv_dir."/".$old[1].".bed")
    {
        $bad_cnv = 1;
    }
    if(-e $bad_snp_bed_dir."/".$old[1].".bed")
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
	
	if ($cancer_type eq "SKCM")
	{
	    if (($old[2] eq $a[2]) & ($old[3] eq '06') & ($a[3] eq '10') | ($a[3] eq '11')) # have matched N
	    {
		unless(($old[3] eq '10') | ($old[3] eq '11') | ($old[3] eq '01'))
		{
		    print PTNO $old[0], "\t", $a[0], "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
		}
	    }
	    elsif (($old[2] eq $a[2]) & ($old[3] eq '01') & ( ($a[3] eq '10') | ($a[3] eq '11'))) # have matched N
	    {
		unless(($old[3] eq '10') | ($old[3] eq '11') | ($old[3] eq '06'))
		{
		    print PTNO $old[0], "\t", $a[0], "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
		}
	    }
	    else
	    {
		print PTNO $old[0], "\t", "NaN", "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
	    }
	}
	elsif (($old[2] eq $a[2]) & ($old[3] eq '01') & ( ($a[3] eq '10') | ($a[3] eq '11'))) # have matched N
        {
            unless (($old[3] eq '10') | ($old[3] eq '11') | ($old[3] eq '06'))
            {
                print PTNO $old[0], "\t", $a[0], "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
            }
        }
        else
        {
            unless (($old[3] eq '10') | ($old[3] eq '11') | ($old[3] eq '06'))
            {
                print PTNO $old[0], "\t", "NaN", "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
            }
        }
	
        @old = @a;
        
       if (-e $bad_cnv_dir."/".$old[1].".bed")
        {
            $bad_cnv = 1;
        }
        else
        {
            $bad_cnv = 0;
        }
        if (-e $bad_snp_bed_dir."/".$old[1].".bed")
        {
            $bad_snp = 1;
        }
        else
        {
            $bad_snp = 0;
        }
    }
    
    if (($old[3] eq '10') | ($old[3] eq '11') | ($old[3] eq '06') & ($cancer_type ne "SKCM"))
    {
	EXIT_IF:
	{
	    last EXIT_IF;
	}
    }
    
    if (($old[3] eq '06') & ($cancer_type eq "SKCM") & (($a[3] eq '10') | ($a[3] eq '11')))
    {
	print PTNO $old[0], "\t", $a[0], "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
    }
    elsif (($old[3] eq '01') & ($cancer_type eq "SKCM") & (($a[3] eq '10') | ($a[3] eq '11')))
    {
	print PTNO $old[0], "\t", $a[0], "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
    }
    elsif (($old[2] eq $a[2]) & ($old[3] eq '01') & (($a[3] eq '10') | ($a[3] eq '11'))) # have matched N
    {
        print PTNO $old[0], "\t", $a[0], "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
    }
    else
    {
        print PTNO $old[0], "\t", "NaN", "\t", substr($old[1],0,16), "\t", $bad_cnv, "\t", $bad_snp, "\n";
    }
    close (PTNI);
    close (PTNO);
}

sub have_cnv_snp
{
    my $cnv_snps_file = shift; #cnv snps file
    my $self = $cnv_snps_file and $cnv_snps_file = shift if ref $cnv_snps_file;
    my $cnv_snps_sort_file = shift; #cnv snps file to sort
    
    open (HCSI,"$cnv_snps_file");
    open (HCSO,">$cnv_snps_sort_file");
    
    while(my $r = <HCSI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        if( ($a[3] == 1) && ($a[4] == 1))
        {
            print HCSO $r, "\n";
        }
    }
    close (HCSI);
    close (HCSO);
}

sub mk_files_snps_cnvs
{
    my $cnv_snps_tn_file = shift; #bad snps and cnv tumor/normal file
    my $self = $cnv_snps_tn_file and $cnv_snps_tn_file = shift if ref $cnv_snps_tn_file;
    my ($ase_path,$gene_level,$duplicates) = @_; #path to directory where ase analysis data is stored, directory where gene level ase data is stored, option to include duplicates
    
    open (MKFI,"$cnv_snps_tn_file");
    open (N,">normal_ff");
    open (T,">tumor_ff");
    
    while(my $r = <MKFI>)
    {
        chomp($r);
        my @a = split("\t",$r);
	if (lc $duplicates eq "n" or lc $duplicates eq "no")
	{
	    $a[2] =~ s/-[0-9]+[a-zA-Z]$//;
	}
	
        unless($a[1] =~ /NaN/i) #prints if there is a matching normal
        {
            print N "$ase_path/$gene_level/".$a[1], "\t", $a[2], "\n";
        }
	print T "$ase_path/$gene_level/".$a[0], "\t", $a[2], "\n";
	
    }
    close (MKFI);
    close (N);
    close (T);
}

sub mk_files_tum_norm
{
    my $ase_file = shift; #file with list of ase files
    my $self = $ase_file and $ase_file = shift if ref $ase_file;
    my ($gene_level_path,$duplicates) = @_; #path to directory where gene level analysis data is stored, option to allow duplicates
    
    open (MKFI,"$ase_file");
    open (N,">normal_ff");
    open (T,">tumor_ff");
    
    while(my $r = <MKFI>)
    {
        chomp($r);
        my @a = split("\t",$r);
	#$a[1] =~ s/-[0-9]+[a-zA-Z]\.[a-zA-Z]+$//;
        
        if($a[4] < 10)#prints if sample is normal
        {
	    if (lc $duplicates eq "y" || lc $duplicates eq "yes")
	    {
		print T "$gene_level_path/".$a[2], "\t", $a[1], "\n";
	    }
	    else
	    {
		print T "$gene_level_path/".$a[2], "\t", $a[3], "\n";
	    }
        }
	else
	{
	    if (lc $duplicates eq "y" || lc $duplicates eq "yes")
	    {
		print N "$gene_level_path/".$a[2], "\t", $a[1], "\n" unless $a[2] =~ /^NaN/i;
	    }
	    else
	    {
		print N "$gene_level_path/".$a[2], "\t", $a[3], "\n" unless $a[2] =~ /^NaN/i;
	    }
	}
    }
    close (MKFI);
    close (N);
    close (T);
}

sub run_problem_cnvs
{
    my $bad_cnv_dir = shift; #directory where bad cnv files are stored
    my $self = $bad_cnv_dir and $bad_cnv_dir = shift if ref $bad_cnv_dir;
    my ($cnv_tn_file,$ase_path,$matrix,$problem_cnvs,$overlapSelect) = @_; #bad snps and cnv tumor/normal file,path to directory where as data is stored, directory that will store matrix data for MatLab etc., directory that will store problem cnv files, overlapSelect command
    
    open (RPCI,"$cnv_tn_file");
    
    mkdir "$ase_path/$problem_cnvs" unless(-d "$ase_path/$problem_cnvs");
    system("rm -f $ase_path/$problem_cnvs/*");
    while(my $r = <RPCI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        
        print "Working on $a[0]\n";
        
        my $bed = "$a[2].bed";
        
        `$overlapSelect $bad_cnv_dir/$bed $ase_path/$matrix/rowlabels.bed $ase_path/cnv_t`;
        pull_problems("$ase_path/cnv_t","$ase_path/$problem_cnvs/$a[0]");
    }
    close (RPCI);
}

sub run_problem_snps
{
    my $bad_snp_dir = shift; #directory where bad snp files are stored
    my $self = $bad_snp_dir and $bad_snp_dir = shift if ref $bad_snp_dir;
    my ($snps_tn_file,$ase_path,$matrix,$problem_snps,$overlapSelect) = @_; #bad snps and cnv tumor/normal file,path to directory where as data is stored, directory that will store matrix data for MatLab etc., directory that will store problem cnv files, overlapSelect command
    
    $bad_snp_dir =~ s/\/$//;
    
    mkdir "$ase_path/$problem_snps" unless(-d "$ase_path/$problem_snps");
    system("rm -f $ase_path/$problem_snps/*");
    
    open (PPSI,"$snps_tn_file");
    
    while(my $r = <PPSI>)
    {
	chomp($r);
	my @a = split("\t",$r);
	
	print STDERR "Working on $a[0]\n";
	
	my $bed = "$a[2].bed";
	
	`$overlapSelect $bad_snp_dir/$bed $ase_path/$matrix/rowlabels.bed $ase_path/snp_t`;
	pull_problems("$ase_path/snp_t","$ase_path/$problem_snps/$a[0]");
    }
    close (PPSI);
}

sub pull_problems
{
    my ($snps_cnvs_pull_file,$snps_cnvs_outdir) = @_; #snps and cnvs file to pull problems, snps and cnvs output directory
    
    open (PPI,"$snps_cnvs_pull_file");
    open (PPO,">$snps_cnvs_outdir");
    while(my $r = <PPI>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print PPO $a[3], "\t", 1, "\n";
    }
    close (PPI);
    close (PPO);
}

#########################################compile_gene_faster subs#############################################
sub Make4colBed
{
   my $sorted_bed_file = shift; #file from the sorted bed directory
   open (my $F, "$sorted_bed_file");
   open (my $out,">$sorted_bed_file.4col");
   while (my $l = <$F>)
   {
        chomp($l);
        my @as = split("\t",$l);
        print STDERR "The bed file $sorted_bed_file contains less than 6 cols!\n" and exit if @as<4;
        print $out join("\t",@as[0..3]),"\n";          
   }
   close $F;
   close $out;
   unlink $sorted_bed_file;
   rename "$sorted_bed_file.4col","$sorted_bed_file",;
}

sub convert2mm
{
    #summarize reads of all four bases (A,C,G,T)and prepare data for TallyOVerDefBed
    my ($ase_count_file,$out_file) = @_;
    
    my $FH = "FH".rand();
    my $OUT = "OUT".rand();
    open ($FH, "$ase_count_file");
    open ($OUT, ">$out_file");
    
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
    my ($snp,$def,$out,$ase_path,$overlapSelect) = @_;
    
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
    open ($OUT_FH,">$out");

    unless ($phred =~ /^[0-9\.]+$/)
    {
        print STDERR "-phred must be a number, not $phred!\n";
        exit;
    }
    unless ($minC =~ /^[0-9]+$/)
    {
        print STDERR "-minC must be a number, not $phred!\n";
        exit;
    }
    
    if($sr == 1)
    {
        #split file
        my $PFH = "P".rand();
        my $MFH = "M".rand();
        my $FFH = "F".rand();
       
        open ($PFH,">$ase_path/plus.$rr.bed");
        open ($MFH,">$ase_path/minus.$rr.bed");
        open ($FFH,$snp);
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
        `$overlapSelect $def $ase_path/plus.$rr.bed -mergeOutput $ase_path/tmp.$rr.plus.out`;
        `$overlapSelect $def $ase_path/minus.$rr.bed -mergeOutput $ase_path/tmp.$rr.minus.out`;
        my $F_FH = "F_FH".rand();
        open ($F_FH,"$ase_path/tmp.$rr.plus.out");
        
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
        open ($FH,"$ase_path/tmp.$rr.minus.out");
        
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
       	`$overlapSelect $def $snp -mergeOutput $ase_path/tmp.$rr.all.out`;
        my $FH_= "FH_".rand();
        open ($FH_,"$ase_path/tmp.$rr.all.out");
     
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
    my ($talliedBed,$RefBed,$gene_ase_out,$ase_tmp_path,$ase_path,$overlapSelect) = @_;
    
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
    
    `$overlapSelect $RefBed $talliedBed $sr $strand -mergeOutput $ase_tmp_path/temp.$rr.bed`;

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
    `sort -k $ff,$ff $ase_tmp_path/temp.$rr.bed > $ase_tmp_path/sorted.$rr.temp.bed`;
    `mv $ase_tmp_path/sorted.$rr.temp.bed $ase_tmp_path/temp.$rr.bed`;
    
    my $BEDFH = "BEDFH".rand();
    open ($BEDFH,"$ase_tmp_path/temp.$rr.bed");
    my $OUTFH = "OUTFH".rand();
    open ($OUTFH,">$out");
    
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
    $tt2 = $x[3];
    $tt2 =~ s/\|/,/g; @tt1 = split(",",$tt2);
    my $print_cd_out = $tt1[0]."|".$tt1[1];#allele1|allele2
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
	    $tt2 =~ s/\|/,/g;
	    @tt1 = split(",",$tt2);
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
            my $ppp = get_p($can,$var,"$ase_path");
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
            
	    #To process data of new transcript id;
            $x[3] =~ s/\|/,/g;
            @counts = split(",",$x[3]);
            $can = $counts[2];
            $var = $counts[3];
            $snps_out = $counts[2]."|".$counts[3];
            $snp_ids_out = $x[0]."-".$x[1];
            $tt2 = $x[3];
	    $tt2 =~ s/\|/,/g;
	    @tt1 = split(",",$tt2);
	    $print_cd_out = $tt1[0]."|".$tt1[1];
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
    my $ppp = get_p($can,$var,"$ase_path");
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
   
    close ($OUTFH);
    close ($BEDFH);
    `rm $ase_tmp_path/temp.$rr.bed`;
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
	    $out = `R --no-save $s $total 0.5 < $RNA_ASE/pbinom.R --quiet|tail -2|head   -1`;
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
	    $out = `R --no-save $t $total 0.5 < $RNA_ASE/pbinom.R --quiet|tail -2|head -1`;
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
    my $p_out = get_p($1,$2,"$RNA_ASE");
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
	    $p_out = get_p($ss,$tt,"$RNA_ASE");
	    if($ss < $tt)
            {
		$p_out = $p_out*-1;
	    }
	}
    }    
    return($out."\t".$p_out);
}
#########################################end of compile_gene_faster subs#############################################

sub ase_filter_overlap
{
    my $self = shift;
    my ($filter_names,$RNA_table_file,$output_names_file) = @_;

    open(RMP,">lookup_file.txt");
    foreach my $rnamp (@$filter_names)
    {
        chomp($rnamp);
        my @split_name = split("\\.",$rnamp);
        my $tcga = $split_name[1];
        $tcga =~ s/-[0-9]+[a-zA-Z]$//; #get rid of the sample id at the end of the TCGA ID
        print RMP $split_name[0],"\t",$split_name[1],"\t",$tcga,"\n";
    }
    close (RMP);

    $parsing->vlookup("lookup_file.txt",3,"$RNA_table_file",2,2,"y","lookup_NaN.txt"); #new $RNA_table file contains only overlap IDs and matches those IDs to the IDs in the rna_look.txt file

    my @filt_mpileups = `cat lookup_NaN.txt`;
    open(MPO,"> $output_names_file");
    foreach my $mpm (@filt_mpileups)
    {
        chomp($mpm);
        if ($mpm !~ /\tNaN+$/i)
        {
            my @new_mp = split("\t",$mpm);
            my $joined_mp = join(".",$new_mp[0],$new_mp[1]); #re-joins the file name after being split and used in vlookup
            print MPO "$joined_mp","\n"; #prints the path to the mpileup file
        }
    }
    close(MPO);      
}

sub print_ase_gene_level_path
{
    my $self = shift;
    my ($ase,$list_dir,$lookup_file,$output_file,$count_gene_archive,$overlap,$RNA_CDS_Mpileups_overlap) = @_;
    
    if (lc $overlap eq "y" or lc $overlap eq "yes")
    {
        `ls $list_dir > $lookup_file`;
        $parsing->vlookup("$ase/$lookup_file",1,"$RNA_CDS_Mpileups_overlap",1,1,"y","$ase/counts_gene_NaN.txt");
        `grep -v NaN $ase/counts_gene_NaN.txt > $ase/counts_gene_pull.txt`;
        $parsing->pull_column("$ase/counts_gene_pull.txt","1","$output_file");
    }
    else
    {
        `ls $ase/$list_dir > $output_file`; 
    }
    
    my @count_gene_files = `cat $ase/$output_file`;
    open (CGA,">$ase/$count_gene_archive");
    foreach my $asea (@count_gene_files)
    {
        print CGA "$ase\t$list_dir/$asea";
    }
    close(CGA);
}

1;
