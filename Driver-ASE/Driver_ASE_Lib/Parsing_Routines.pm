package Driver_ASE_Lib::Parsing_Routines;

use strict;
use warnings;
use File::Find;
use File::Copy;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

no strict "refs";

###############################################################################
# Parsing_Routines.pm version 1
# Written by Tomas Babak, Zhongshan Cheng, Michael Vermeulen and Micheal Rollins-Green
#  in the Department of Biology at Queen's University.
#
# Parsing_Routines
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

###############################################################################
# Program Usage
###############################################################################
sub usage
{
    my $script = shift;
    my $self = $script and $script = shift if ref $script;
    
    if(!defined $script)
    {
        $script = "";
    }
    
    if($script eq "0")
    {
        print "usage: program [--disease_abbr|-d disease_abbr (e.g. PRAD)] [--exp_strat|-e Experimental Strategy (Genotyping array)] [--array_type|-a array data type (e.g. Genotypes)] [--command|-c curl or aria2c] [--key|-k path to gdc key] [--help|-h]\n";
        
        print"Any names with spaces must be wraped in DOUBLE QUOTES or have back slashes to escape spaces.\n";
        
        print 'Experimental Strategies used:
        Copy number estimate
        Genotypes', "\n";
    }
    if($script eq "0_table")
    {
        print "usage: program [--disease_abbr|-d disease_abbr (e.g. PRAD)] [--exp_strat|-e Experimental Strategy (Genotyping array)] [--array_type|-a array data type (e.g. Genotypes)] [--help|-h]\n";
        
        print"Any names with spaces must be wraped in DOUBLE QUOTES or have back slashes to escape spaces.\n";
        
        print 'Experimental Strategies used:
        Copy number estimate
        Genotypes', "\n";
    }
    elsif($script eq "3.0")
    {
        print "usage: program [--disease_abbr|-d disease_abbr (e.g. PRAD)] [--exp_strat|-e Experimental Strategy (e.g. WGS/RNA-Seq)] [--option|-o all, download or mpileups] [--number|-n number of bams to download for RNA-Seq and number of bam pairs for WGS] [--command|-c curl or aria2c] [--var_path|-v path to the VarScan jar file (Enter if -e is WGS)] [--key|-k peht to gdc key] [--help|-h]\n";
    }
    elsif($script eq "3.0_table")
    {
        print "usage: program [--disease_abbr|-d disease_abbr (e.g. PRAD)] [--exp_strat|-e Experimental Strategy (e.g. WGS/RNA-Seq)] [--command|-c curl or aria2c] [--key|-k peht to gdc key] [--help|-h]\n";
    }
    elsif($script eq "3.2")
    {
	print "usage: script [--disease_abbr|-d disease_abbr (e.g. PRAD)] [--bsc|-b (y|yes|n|no)] [--overlap|-o (y|yes|n|no)] [--duplicates|-du (y|yes|n|no)] [--compress|-c (y|yes|n|no)] [--remfiles|-r (y|yes|n|no)] [--help|-h]\n";
    }
    elsif($script eq "4.0")
    {
	print "usage: program [--disease_abbr|-d (e.g. PRAD)] [--readcutoff|-r read cutoff (e.g. 20)] [--tfreq|-t tumor alt frequency (e.g. 0.1)] [--nfreq|-n normal alt frequency (e.g. 0.1)]"
    }
    else
    {
        print "usage: script [--disease_abbr|-d disease_abbr (e.g. PRAD)] [--help|-h]\n";
    }
    
    print '
        All available TCGA cancer types:
	
	Breast invasive carcinoma [BRCA]
        Glioblastoma multiforme [GBM]
        Ovarian serous cystadenocarcinoma [OV]
        Lung adenocarcinoma [LUAD]
        Uterine Corpus Endometrial Carcinoma [UCEC]
        Kidney renal clear cell carcinoma [KIRC]
        Head and Neck squamous cell carcinoma [HNSC]
        Brain Lower Grade Glioma [LGG]
        Thyroid carcinoma [THCA]
        Lung squamous cell carcinoma [LUSC]
        Prostate adenocarcinoma [PRAD]
        Stomach adenocarcinoma [STAD]
        Skin Cutaneous Melanoma [SKCM]
        Colon adenocarcinoma [COAD]
        Bladder Urothelial Carcinoma [BLCA]
        Liver hepatocellular carcinoma [LIHC]
        Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC]
        Kidney renal papillary cell carcinoma [KIRP]
        Sarcoma [SARC]
        Esophageal Carcinoma [ESCA]
        Pancreatic Adenocarcinoma [PAAD]
        Pheochromocytoma and Paraganglioma [PCPG]
        Rectum Adenocarcinoma [READ]
        Testicular Germ Cell Tumors [TGCT]
        Acute Myeloid Leukemia [LAML]
        Thymoma [THYM]
        Kidney Chromophobe [KICH]
        Adernocortical Carcinoma [ACC]
        Mesothelioma [MESC]
        Uveal Melanoma [UVM]
        Lymphoid Neoplasm Diffuse Large B-cell Lymphoma [DLBC]
        Uterine Carcinosarcoma [UCS]
        Cholangiocarcinoma [CHOL]',"\n";
    exit;
}

###############################subs#############################

sub copyfile2newfullpath
{
    my $gdc_key_fullpath = shift;
    my $self = $gdc_key_fullpath and $gdc_key_fullpath = shift if ref $gdc_key_fullpath;
    my $newfullpath = shift;
    #check the existense of new dir of the newfullpath;
    my ($newdir) = "$newfullpath" =~ /^(.*)\/[^\/]+/;
    `mkdir -p "$newdir"` unless (-d $newdir);
    if(-f $gdc_key_fullpath)
    {   
        `cp -f $gdc_key_fullpath \'$newfullpath\'`;
        print "Now copying the key: $gdc_key_fullpath to the new fullpath $newfullpath $!\n";
    }
    else
    {
       print STDERR "The path to the gdc key doesn\'t exist:\n","$gdc_key_fullpath\n";
       exit;      
    }
}

sub get_only_files_in_dir
{
    my $path = shift;
    my $self = $path and $path = shift if ref $path;
    opendir my $dir, "$path" or die("Can't open the directory $path: $!\n");
    my @get_files = readdir $dir;
    close($dir);    

    @get_files = grep{-f "$path/$_"}@get_files; 
    
    return @get_files;
}

sub QueryBase
{
    my $Qfile = shift;
    my $self = $Qfile and $Qfile = shift if ref $Qfile;
    my $queryCol = shift;
    my $sourceFile = shift;
    my $returnCols_rgx = shift;
    my $outfile = shift;
    my $KeepAllColsinQeryFile = shift;
    
    $queryCol = $queryCol - 1;
    
    open(BASE,$sourceFile) or die "Cannot open base file: $sourceFile: $!\n";
    open(Q,$Qfile) or die "Cannot open query file: $Qfile: $!\n";
    open (RES,">$outfile") or die "Can't open output file: $outfile for output: $!\n";
    my %hash;
    while(my $l = <Q>)
    {
    	chomp($l);
    	my @line1;
    	if($l =~ /^\s*$/)
        {
            next;
	}
        else
        {
	    @line1 = split("\t",$l);
            if($KeepAllColsinQeryFile eq 1)
            {
                push @{$hash{$line1[$queryCol]}},$l unless $line1[$queryCol] eq "";
            }
            else
            {
                push @{$hash{$line1[$queryCol]}},$line1[$queryCol] unless $line1[$queryCol] eq "";   
            }
        }
    } 
    close Q;
    
    chomp(my @temp_array = <BASE>);
    close(BASE);
    my @headers = split("\t",$temp_array[0]);
    foreach my $q(keys %hash)
    {
        foreach my $b (@temp_array)
        {
	    my @tmp = split("\t",$b);
            my $match_idx = get_element_position_in_array($q,\@tmp);
            my $rgx_match_idx = get_rgx_matched_position_in_array($returnCols_rgx,\@tmp);
            if($match_idx >= 0 and $rgx_match_idx >= 0)
            {
		foreach my $v (@{$hash{$q}})
                {
                    print RES $v,"\t",$tmp[$rgx_match_idx],"\n";
		}
            }
            elsif($match_idx >= 0)
            {
		foreach my $v (@{$hash{$q}})
                {
                    print RES $v,"\t","NaN","\n";
		}
	    }
        }
    }
   close(RES);
 }

sub get_element_position_in_array
{
    my $search_for = shift;
    my $self = $search_for and $search_for = shift if ref $search_for;
    my $target_array_ref = shift;
    my ($index) = grep {$$target_array_ref[$_] eq $search_for}0..$#$target_array_ref;
    
    if(defined($index))
    {
        return $index;
    }
    else
    {
        return -500;
    }
}

sub get_rgx_matched_position_in_array
{
    my $search_for = shift;
    my $self = $search_for and $search_for= shift if ref $search_for;
    my $target_array_ref = shift;
    my ($index) = grep { $$target_array_ref[$_] =~ /$search_for/i } 0..$#$target_array_ref;
    
    if(defined($index))
    {
        return $index;
    }
    else
    {
        return -500;
    }
}

sub vlookup
{
    my $lookupfile = shift;
    my $self = $lookupfile and $lookupfile = shift if ref $lookupfile;
    my $queryCol = shift;
    my $sourceFile = shift;
    my $lookupCol = shift;
    my $returnCols = shift;
    my $append = shift;
    my $outfile = shift;
    
    open(L,$sourceFile) or die "Cannot open $sourceFile: $!\n";
    open (VO,">$outfile") or die "Can't open $outfile for output: $!\n";
    # load required columns into memory
    my @cols = split(",",$returnCols);
    my $col = $lookupCol;
    $col -= 1;
    my $nan = '';
    for(my $i = 0;$i < scalar(@cols);$i++)
    {
        $cols[$i] -= 1; $nan = $nan."\t"."NaN"
    }
    $nan =~ s/^\t//;
    
    my %look = ();
    while(my $r = <L>)
    {
        chomp($r);
        if($r eq "")
        {
            next;
        }
        else
        {
            my @a = split("\t",$r);
            $look{$a[$col]} = join("\t",@a[@cols]);    
        }
    }
    close(L);
    
    # parse through sourceFile and print, appending if $append eq y
    open(S,$lookupfile) or die "Cannot open $lookupfile: $!\n";
    while(my $r = <S>)
    {
        chomp($r);
        if($r eq '')
        {
            next;
        }
        my @a = split("\t",$r);
        if(exists($look{$a[$queryCol-1]}))
        { 
           if($append eq 'y')
            {
                print VO $r, "\t", $look{$a[$queryCol-1]}, "\n";
            }
           else
            {
                print VO $look{$a[$queryCol-1]}, "\n";
            }
        }
        else
        {
            #not in sourceFile
            if($append eq 'y')
            {
                print VO $r, "\t", $nan, "\n";
            }
            else
            {
                print VO $nan, "\n";
            }
        }
    }
    close(VO);
    close(S);
}

sub pull_column
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my @cols = split(",",shift);
    my $outfile = shift;
    
    open(P,$infile) or die("Can't open $infile for input. $!\n");
    open(C,">$outfile") or die("Can't open $outfile for output. $!\n");
    for(my $i=0;$i < scalar(@cols);$i++)
    {
        $cols[$i]-=1;
    }
    
    while(my $r = <P>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        print C join("\t",@a[@cols]), "\n";
    }
    close(P);
    close(C);
}

sub strip_head
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    my $num_lines = shift || 1;
    my $col_nums = shift;
    my $sep = shift;
    
    open(S,$infile) or die("Can't open file $infile for input: $!\n");
    open(H,">$outfile") or die("Can't open file $outfile for output");
    for(my $i = 0;$i < $num_lines;$i++)
    {
        my $r = <S>;
    }
    while(my $r = <S>)
    {
        if(defined $col_nums)
        {
            my @tmps;
            chomp($r);
            @tmps = split("\t",$r);
            @tmps = split("$sep",$r) if defined $sep;
            if(defined $sep)
            {
                print H join("$sep",$tmps[$col_nums]),"\n"; 
            }
            else
            {
                print H join("\t",$tmps[$col_nums]),"\n"; 
            }  
        }
        else
        {            
            print H $r;            
        }
    }
    close(S);
    close(H);
}

sub matricize
{
    my $index_file = shift;
    my $self = $index_file and $index_file = shift if ref $index_file;
    my $print_file = shift;
    my $idx_col_print = shift;
    my $out_col = shift;
    my $dirpath = shift;#wgs path or rna path
    
    # get all ids and unique them
    print STDERR "getting identifiers\n";
    open(F,$index_file);
    my @ids;
    
    while (my $r = <F>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        open(FF,$a[0]) or die "$a[0] does not exist\n";
        while (my $rr=<FF>)
        {
            chomp($rr);
            my @aa = split("\t",$rr);
            unshift(@ids,$aa[$idx_col_print-1]);
        }
        close(FF);
    }
    close(F);
    @ids=unique(@ids);
    
    # for each file, get it into a hash, then append to an @out array
    print STDERR "loading in data\n";
    my @out=();
    my @samples;
    my $j=0;
    open(F,$print_file);
    while (my $r = <F>) 
    {
        chomp($r);
        my @a = split("\t",$r);
            push(@samples,$a[1]);
        print STDERR "working on $a[1]\n";
            open(FF,$a[0]) or die "$a[1] does not exist\n";
            my %hash=();
        while (my $rr=<FF>)
        {
            chomp($rr);
            my @aa = split("\t",$rr);
            $hash{$aa[$idx_col_print-1]}=$aa[$out_col-1];
        }
        close(FF);
       
        # append to out
        for(my $i=0; $i<scalar(@ids); ++$i)
        {
            if(exists($hash{$ids[$i]}))
            {
                $out[$j][$i]=$hash{$ids[$i]};
            }
            else
            {
                $out[$j][$i]='NaN';
            }
        }
        ++$j;
       
    }
    close(F);
    
    # print contents out
    print STDERR "printing output\n";
    open(COL,">$dirpath/collabels.txt");
    print COL join("\n",@samples), "\n";
    close(COL);
    open(ROW,">$dirpath/rowlabels.txt");
    print ROW join("\n",@ids), "\n";
    close(ROW);
    open(MAT,">$dirpath/matrix.tab");
    for(my $i=0; $i<scalar(@ids); ++$i)
    {
        my $lineout='';
        for(my $k=0; $k<$j; ++$k)
        {
            $lineout .= $out[$k][$i]."\t";
        }
        $lineout=~s/\t$/\n/;
        print MAT $lineout;
    }
    close(MAT);
}

sub matricize_version_two
{
    my $index_file = shift;
    my $self = $index_file and $index_file = shift if ref $index_file;
    my $print_file = shift;
    my $index_column = shift;
    my $print_index_column = shift;
    my $print_output_column = shift;
    
    # get all ids and unique them
    print STDERR "getting identifiers\n";
    open(F,$index_file);
    my @ids;
    
    while (my $r = <F>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        open(FF,$a[0]) or die "$a[0] does not exist\n";
        while (my $rr = <FF>)
        {
            chomp($rr);
            my @aa = split("\t",$rr);
            unshift(@ids,$aa[$index_column-1]);
        }
        close(FF);
    }
    close(F);
    @ids = unique(@ids);
    
    # for each file, get it into a hash, then append to an @out array
    print STDERR "loading in data\n";
    my @out = ();
    my @samples;
    my $j = 0;
    open(F,$print_file);
    while (my $r = <F>) 
    {
        chomp($r);
        my @a = split("\t",$r);
            push(@samples,$a[1]);
        print STDERR "working on $a[1]\n";
            open(FF,$a[0]) or die "$a[1] does not exist\n";
            my %hash = ();
        while (my $rr=<FF>)
        {
            chomp($rr);
            my @aa = split("\t",$rr);
            $hash{$aa[$print_index_column-1]} = $aa[$print_output_column-1];
        }
        close(FF);
        
        # append to out
        for(my $i = 0;$i < scalar(@ids);$i++)
        {
            if(exists($hash{$ids[$i]}))
            {
                $out[$j][$i] = $hash{$ids[$i]};
            }
            else
            {
                $out[$j][$i] = 'NaN';
            }
        }
        ++$j;   
    }
    close(F);
    
    # print contents out
    print STDERR "printing output\n";
    open(OUT,">collabels.txt");
    print OUT join("\n",@samples), "\n";
    close(OUT);
    open(OUT,">rowlabels.txt");
    print OUT join("\n",@ids), "\n";
    close(OUT);
    open(OUT,">matrix.tab");
    for(my $i = 0;$i < scalar(@ids);$i++)
    {
        my $lineout = '';
        for(my $k = 0;$k < $j;$k++)
        {
            $lineout .= $out[$k][$i]."\t";
        }
        $lineout =~ s/\t$/\n/;
        print OUT $lineout;
    }
    close OUT;
}

############matricize and matricize.v2 sub####################
sub unique 
{
    my (@in) = @_;
    my %unique = ();
    foreach (@in)
    {
        $unique{$_} = 1;
    }
    return(keys %unique);
}
############end of matricize and matricize_version_two sub####################

sub Get_PIDs
{
    my $self = shift;
    my $exclude_jobs_ref;
    if (ref $self)
    {
        $exclude_jobs_ref = \@_;
    }
    else
    {
        $exclude_jobs_ref = $self;
    }
    
    if (scalar(@$exclude_jobs_ref) < 1)
    {
       $exclude_jobs_ref = [0];
    }
   
    my $FH=FileHandle->new("ps a |");
    my @jobs;
    while (my $pid=<$FH>)
    {
       chomp($pid);
       $pid=~s/^\s+//g;
       if ($pid =~ /^\d+/)
        {
            my @as=split("\t| ",$pid);
            my $job = $as[0];
            #print $job,"\n";
            my $t = 0;
            foreach my $j(@$exclude_jobs_ref)
            {
                $t++ if ($job == $j);
            }
            push @jobs,$job if $t == 0;
        }
    }
    $FH->close;
    print join(" ",@jobs),"\n" if @jobs > 0;
}

sub rand_num
{
    my $rand = rand();
    my $rd_num = substr($rand,2,6);
    return($rd_num);   
}

sub temp_filename
{
    my $tmp = "tmp".&rand_num;
    return($tmp);
}

sub get_alt
{
    my ($self,$in) = @_;

    my @aa;

    $aa[0] = ($in =~ tr/aA/aA/);
    $aa[1] = ($in =~ tr/cC/cC/);
    $aa[2] = ($in =~ tr/gG/gG/);
    $aa[3] = ($in =~ tr/tT/tT/);

    my @list_order = sort{$aa[$b] <=> $aa[$a]}0..$#aa;
    my $acgt_sum = $aa[0] + $aa[1] + $aa[2] +$aa[3];
    
    my $prop = 0;
    if($acgt_sum > 0)
    {
        $prop = $aa[$list_order[1]]/$acgt_sum;
    }
    return(substr('ACGT',$list_order[1],1),$prop);
}

sub get_ref
{
    my ($self,$in) = @_;

    my @aa;

    $aa[0] = ($in =~ tr/aA/aA/);
    $aa[1] = ($in =~ tr/cC/cC/);
    $aa[2] = ($in =~ tr/gG/gG/);
    $aa[3] = ($in =~ tr/tT/tT/);

    my @list_order = sort {$aa[$b] <=> $aa[$a]}0..$#aa;
    my $acgt_sum = $aa[0] + $aa[1] + $aa[2] + $aa[3];
    
    my $prop = 0;
    if($acgt_sum > 0)
    {
        $prop = $aa[$list_order[1]]/$acgt_sum;
    }
    return(substr('ACGT',$list_order[0],1),$prop);
}

1;
