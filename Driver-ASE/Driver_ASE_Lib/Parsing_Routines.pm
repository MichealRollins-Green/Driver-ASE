package Driver_ASE_Lib::Parsing_Routines;

use strict;
use warnings;
use File::Find;
use File::Copy;
use File::Which;
use Archive::Tar;
use autodie;
no warnings 'once';

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
    
    if (!defined $script)
    {
        $script = "";
    }
    
    if ($script eq "0")
    {
        print "usage: program [--cancer|-c cancer_type (e.g. PRAD)] [--Expstrategy|-E Experimental strategy (Genotyping array)] [--arraytype|-a array data type (e.g. Genotypes or “Copy number estimate”)] [--download|-d curl or aria2c] [--key|-k path to gdc key] [--help|-h]\n";
        
        print"Any names with spaces must be wraped in DOUBLE QUOTES or have back slashes to escape spaces.\n";
        
        print 'Experimental Strategies used:
        Copy number estimate
        Genotypes', "\n";
    }
    elsif ($script eq "0_table")
    {
        print "usage: program [--cancer|-c cancer_type (e.g. PRAD)] [--Expstrategy|-E Experimental strategy (Genotyping array)] [--arraytype|-a array data type (e.g. Genotypes or “Copy number estimate”)] [--help|-h]\n";
        
        print"Any names with spaces must be wraped in DOUBLE QUOTES or have back slashes to escape spaces.\n";
        
        print 'Experimental Strategies used:
        Copy number estimate
        Genotypes', "\n";
    }
    elsif ($script eq "1.1")
    {
	print "usage: script [--cancer|-c cancer_type (e.g. PRAD)] [--plink|-p path to plink] [--sampletype|-s 0(normal)|1(tumor)|2(normal/tumor)] [--help|-h]\n";
    }
    elsif ($script eq "1.2")
    {
	print "usage: script [--cancer|-c cancer_type (e.g. PRAD)] [--shapeit|-s path to shapeit] [--impute|-i path to impute2] [--help|-h]\n";
    }
    elsif ($script eq "1.3")
    {
	print "usage: script [--cancer|-c cancer_type (e.g. PRAD)] [--plink|-p path to plink] [--help|-h]\n";
    }
    elsif ($script eq "3.0")
    {
        print "usage: program [--cancer|-c cancer_type (e.g. PRAD)] [--Expstrategy|-E Experimental strategy (e.g. WGS/RNA-Seq)] [--option|-o all, download or mpileups] [--intersect|-i (y|yes|n|no)] [--number|-n number of bams to download for RNA-Seq and number of bam pairs for WGS] [--download|-d curl or aria2c] [--samtools|-s path to samtools] [--VarScan|-V path to the VarScan jar file (Enter if -E is WGS)] [--key|-k path to gdc key] [--help|-h]\n";
    }
    elsif ($script eq "3.0_table")
    {
        print "usage: program [--cancer|-c cancer_type (e.g. PRAD)] [--Expstrategy|-E Experimental strategy (e.g. WGS/RNA-Seq)] [--download|-d curl or aria2c] [--overlap|-o (n|no|y|yes)] [--key|-k path to gdc key] [--help|-h]\n";
    }
    elsif ($script eq "3.1")
    {
	print "usage: script [--cancer|-c cancer_type (e.g. PRAD)] [--overlapSelect|-S path to overlapSelect] [--overlap|-o (y|yes|n|no)] [--help|-h]\n";
    }
    elsif ($script eq "3.2")
    {
	print "usage: script [--cancer|-c cancer_type (e.g. PRAD)] [--overlapSelect|-S path to overlapSelect] [--badsnpcnvs|-b (y|yes|n|no)] [--overlap|-o (y|yes|n|no)] [--duplicates|-d (y|yes|n|no)] [--remfiles|-r (y|yes|n|no)] [--archive|-a (y|yes|n|no)] [--help|-h]\n";
    }
    elsif ($script eq "4.0")
    {
	print "usage: program [--cancer|-c (e.g. PRAD)] [--readcutoff|-r read cutoff (e.g. 20)] [--tfreq|-t tumor alt frequency (e.g. 0.1)] [--nfreq|-n normal alt frequency (e.g. 0.1)] [--overlap|-o (y|yes|n|no)] [--help|-h]\n"
    }
    elsif ($script eq "4.1")
    {
	print "usage: script [--cancer|-c cancer_type (e.g. PRAD)] [--overlap|-o (y|yes|n|no)] [--remfiles|-r (y|yes|n|no)] [--archive|-a (y|yes|n|no)] [--overlapSelect|-S path to overlapSelect] [--help|-h]\n"
    }
    else
    {
        print "usage: script [--cancer|-c cancer_type (e.g. PRAD)] [--help|-h]\n";
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

sub check_directory_existence
{
    my $self = shift; #assigns reference
    my @directories = @_; #directories and/or files to search for
    my @msg;
    my $cnt = 0;
    foreach my $dir (@directories)
    {
	if (!(-d $dir) and !(-f $dir))
	{
	    print STDERR "$dir does not exist! It was either deleted, moved or renamed.\n";
	    if ($dir =~ /Database$/i)
	    {
		$msg[$cnt] = "If the missing directory is the Database directory, it may not have been downloaded or the name was changed. Check the README.md file on the github page to find out where to get the Database directory.\n";
		$cnt += 1;
	    }
	    else
	    {
		$msg [$cnt] = "If $dir is missing, It is likely that it was deleted, moved, renamed or the scripts were ran out of order.\n";
		$cnt += 1;
	    }
	}
	else
	{
	    print "Found directory/file $dir.\n";
	}
    }
    
    if (scalar(@msg) > 0)
    {
	foreach my $df (@msg)
	{
	    print STDERR "$df";
	}
	exit;
    }
}

sub determine_option
{
    my $self = shift;
    my ($RNA_Path,$cds_sorted,$option) = @_; #path to directory where RNA-Seq data from analysis is stored, path to sorted cds fles, download option
    if (!(-d "$RNA_Path/$cds_sorted"))
    {
	print "$RNA_Path/$cds_sorted does not exist! Only downloads will be performed!\n";
	return $option = "download";
    }
    else
    {
	print "Found directory $RNA_Path/$cds_sorted.\n";
	return $option;
    }
}

sub check_cancer_type
{
    my $database_path = shift;
    my $self = $database_path and $database_path = shift if ref $database_path;
    my $cancer_type = shift;
    
    #Check if the cancer type entered exists with in the file.
    open (my $can,"$database_path/Cancer_Types.txt");
    my @can_types = <$can>;
    my $line_num = 1;
    my $num_of_ctypes = scalar(@can_types);
    my $no_count = 0;
    
    print "Locating $cancer_type...\n";
    
    foreach my $line (@can_types)
    {
	chomp($line);
	if ($cancer_type eq $line)
	{
	    print "Found $cancer_type on line $line_num.\n\nContinuing program.\n\n";
	    return $cancer_type;
	}
	else
	{
	    print "No $cancer_type on line $line_num.\n";
	    print "Line $line_num was $line.\n\n";
	    $no_count += 1;
	}
	$line_num += 1;
    }
    close ($can);
    
    if ($no_count == $num_of_ctypes)
    {
	print "$cancer_type is not in the $database_path/Cancer_Types.txt file! Maybe it was misspelled or it does not exist within the file.\n";
	print "Skipping to the next possible cancer type entered.\n";
	return undef;
    }
}

sub CheckSoftware
{
    my $self = shift;
    my @software = @_;
    
    for my $soft(@software)
    {
	#check if only the name of the software was entered
	if ($soft =~ /^[a-zA-Z|0-9]+$/)
	{
	    my $which = which $soft;
	    if (!defined $which)
	    {
		print STDERR "$soft does not exist or is not in the global path!\n";
		print STDERR "Make sure that the name of the software is spelled correctory and check where it might be located.\n";
		exit;
	    }
	    else
	    {
		print "$soft found! The full path is: $which\n";
	    }
	}
	else
	{
	    #if software name was not entered, check if the full path to the software was entered
	    if (!(-f "$soft"))
	    {
		print STDERR "$soft does not exist!\n";
		print STDERR "Make sure that the full path has been input correctly.\n";
		exit;
	    }
	    else
	    {
		print "$soft found\n";
	    }
	}
    }
}

sub copyfile2newfullpath
{
    my $gdc_key_location = shift;
    my $self = $gdc_key_location and $gdc_key_location = shift if ref $gdc_key_location;
    my $gdc_copy_location = shift;
    #check the existense of new dir of the newfullpath;
    my ($newdir) = "$gdc_copy_location" =~ /^(.*)\/[^\/]+/;
    `mkdir -p "$newdir"` unless (-d $newdir);
    `cp -f $gdc_key_location \'$gdc_copy_location\'`;
    print "Now copying $gdc_key_location to the new fullpath $gdc_copy_location.\n\n$!\n";
}

sub get_only_files_in_dir
{
    my $path = shift;
    my $self = $path and $path = shift if ref $path;
    opendir (my $dir, "$path");
    my @get_files = readdir $dir;
    closedir ($dir);    

    @get_files = grep{-f "$path/$_"}@get_files; 
    
    return @get_files;
}

sub QueryBase
{
    my $Qfile = shift;
    my $self = $Qfile and $Qfile = shift if ref $Qfile;
    my ($queryCol,$sourceFile,$returnCols_rgx,$outfile,$KeepAllColsinQeryFile) = @_;
    
    $queryCol = $queryCol - 1;
    
    open (BASE,$sourceFile);
    open (Q,$Qfile);
    open (RES,">$outfile");
    my %hash;
    while (my $l = <Q>)
    {
    	chomp($l);
    	my @line1;
    	if ($l =~ /^\s*$/)
        {
            next;
	}
        else
        {
	    @line1 = split("\t",$l);
            if ($KeepAllColsinQeryFile eq 1)
            {
                push @{$hash{$line1[$queryCol]}},$l unless $line1[$queryCol] eq "";
            }
            else
            {
                push @{$hash{$line1[$queryCol]}},$line1[$queryCol] unless $line1[$queryCol] eq "";   
            }
        }
    } 
    close (Q);
    
    chomp(my @temp_array = <BASE>);
    close (BASE);
    my @headers = split("\t",$temp_array[0]);
    foreach my $q (keys %hash)
    {
        foreach my $b (@temp_array)
        {
	    my @tmp = split("\t",$b);
            my $match_idx = get_element_position_in_array($q,\@tmp);
            my $rgx_match_idx = get_rgx_matched_position_in_array($returnCols_rgx,\@tmp);
            if ($match_idx >= 0 and $rgx_match_idx >= 0)
            {
		foreach my $v (@{$hash{$q}})
                {
                    print RES $v,"\t",$tmp[$rgx_match_idx],"\n";
		}
            }
            elsif ($match_idx >= 0)
            {
		foreach my $v (@{$hash{$q}})
                {
                    print RES $v,"\t","NaN","\n";
		}
	    }
        }
    }
   close (RES);
 }

sub get_element_position_in_array
{
    my $search_for = shift;
    my $self = $search_for and $search_for = shift if ref $search_for;
    my $target_array_ref = shift;
    my ($index) = grep {$$target_array_ref[$_] eq $search_for}0..$#$target_array_ref;
    
    if (defined($index))
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
    my $self = $search_for and $search_for = shift if ref $search_for;
    my $target_array_ref = shift;
    my ($index) = grep { $$target_array_ref[$_] =~ /$search_for/i } 0..$#$target_array_ref;
    
    if (defined($index))
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
    my ($queryCol,$sourceFile,$lookupCol,$returnCols,$append,$outfile) = @_;
    
    open (L,$sourceFile);
    open (VO,">$outfile");
    # load required columns into memory
    my @cols = split(",",$returnCols);
    my $col = $lookupCol;
    $col -= 1;
    my $nan = '';
    for (my $i = 0;$i < scalar(@cols);$i++)
    {
        $cols[$i] -= 1; $nan = $nan."\t"."NaN"
    }
    $nan =~ s/^\t//;
    
    my %look = ();
    while (my $r = <L>)
    {
        chomp($r);
        if ($r eq "")
        {
            next;
        }
        else
        {
            my @a = split("\t",$r);
            $look{$a[$col]} = join("\t",@a[@cols]);    
        }
    }
    close (L);
    
    # parse through sourceFile and print, appending if $append eq y
    open (S,$lookupfile);
    while (my $r = <S>)
    {
        chomp($r);
        if ($r eq '')
        {
            next;
        }
        my @a = split("\t",$r);
        if (exists($look{$a[$queryCol-1]}))
        { 
           if ($append eq 'y')
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
            if ($append eq 'y')
            {
                print VO $r, "\t", $nan, "\n";
            }
            else
            {
                print VO $nan, "\n";
            }
        }
    }
    close (VO);
    close (S);
}

sub pull_column
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my @cols = split(",",shift);
    my $outfile = shift;
    
    open (P,$infile);
    open (C,">$outfile");
    for (my $i=0;$i < scalar(@cols);$i++)
    {
        $cols[$i]-=1;
    }
    
    while (my $r = <P>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        print C join("\t",@a[@cols]), "\n";
    }
    close (P);
    close (C);
}

sub strip_head
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    my $num_lines = shift || 1;
    my $col_nums = shift;
    my $sep = shift;
    
    open (S,$infile);
    open (H,">$outfile");
    for (my $i = 0;$i < $num_lines;$i++)
    {
        my $r = <S>;
    }
    while (my $r = <S>)
    {
        if (defined $col_nums)
        {
            my @tmps;
            chomp($r);
            @tmps = split("\t",$r);
            @tmps = split("$sep",$r) if defined $sep;
            if (defined $sep)
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
    close (S);
    close (H);
}

sub matricize
{
    my $index_file = shift;
    my $self = $index_file and $index_file = shift if ref $index_file;
    my ($print_file,$idx_col_print,$out_col,$dirpath) = @_;
    
    # get all ids and unique them
    print "Getting identifiers...\n";
    open (F,$index_file);
    my @ids;
    
    while (my $r = <F>) 
    {
        chomp($r);
        my @a = split("\t",$r);
       open (FF,$a[0]);
        while (my $rr=<FF>)
        {
            chomp($rr);
            my @aa = split("\t",$rr);
            unshift(@ids,$aa[$idx_col_print-1]);
        }
        close (FF);
    }
    close (F);
    @ids=unique(@ids);
    
    # for each file, get it into a hash, then append to an @out array
    print "Loading in data...\n";
    my @out = ();
    my @samples;
    my $j = 0;
    open (F,$print_file);
    while (my $r = <F>) 
    {
        chomp($r);
        my @a = split("\t",$r);
            push(@samples,$a[1]);
        print "Working on $a[1]\n";
           open (FF,$a[0]);
            my %hash = ();
        while (my $rr = <FF>)
        {
            chomp($rr);
            my @aa = split("\t",$rr);
            $hash{$aa[$idx_col_print-1]} = $aa[$out_col-1];
        }
        close (FF);
       
        # append to out
        for (my $i = 0;$i < scalar(@ids);$i++)
        {
            if (exists($hash{$ids[$i]}))
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
    close (F);
    
    # print contents out
    print "Printing output...\n";
    open (COL,">$dirpath/collabels.txt");
    print COL join("\n",@samples), "\n";
    close (COL);
    open (ROW,">$dirpath/rowlabels.txt");
    print ROW join("\n",@ids), "\n";
    close (ROW);
    open (MAT,">$dirpath/matrix.tab");
    for (my $i = 0;$i < scalar(@ids);$i++)
    {
        my $lineout = '';
        for (my $k = 0;$k < $j;$k++)
        {
            $lineout .= $out[$k][$i]."\t";
        }
        $lineout=~s/\t$/\n/;
        print MAT $lineout;
    }
    close (MAT);
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
    print "Getting identifiers...\n";
    open (F,$index_file);
    my @ids;
    
    while (my $r = <F>) 
    {
        chomp($r);
        my @a = split("\t",$r);
	open (FF,$a[0]);
        while (my $rr = <FF>)
        {
            chomp($rr);
            my @aa = split("\t",$rr);
            unshift(@ids,$aa[$index_column-1]);
        }
        close (FF);
    }
    close (F);
    @ids = unique(@ids);
    
    # for each file, get it into a hash, then append to an @out array
    print "Loading in data...\n";
    my @out = ();
    my @samples;
    my $j = 0;
    open (F,$print_file);
    while (my $r = <F>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        push(@samples,$a[1]);
        print "Working on $a[1]\n";
        open (FF,$a[0]);
        my %hash = ();
        while (my $rr = <FF>)
        {
            chomp($rr);
            my @aa = split("\t",$rr);
            $hash{$aa[$print_index_column-1]} = $aa[$print_output_column-1];
        }
        close (FF);
        
        # append to out
        for (my $i = 0;$i < scalar(@ids);$i++)
        {
            if (exists($hash{$ids[$i]}))
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
    close (F);
    
    # print contents out
    print "Printing output...\n";
    open (OUT,">collabels.txt");
    print OUT join("\n",@samples), "\n";
    close (OUT);
    open (OUT,">rowlabels.txt");
    print OUT join("\n",@ids), "\n";
    close (OUT);
    open (OUT,">matrix.tab");
    for (my $i = 0;$i < scalar(@ids);$i++)
    {
        my $lineout = '';
        for (my $k = 0;$k < $j;$k++)
        {
            $lineout .= $out[$k][$i]."\t";
        }
        $lineout =~ s/\t$/\n/;
        print OUT $lineout;
    }
    close (OUT);
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

sub Overlap_RNA_WGS_Geno
{
    my $Analysispath = shift;
    my $self = $Analysispath and $Analysispath = shift if ref $Analysispath;
    my ($RNA_table,$WGS_table,$Geno_table,$RNA_overlap,$WGS_overlap,$Geno_overlap,$Intersect,$Exp_Strategy,$cancer_type) = @_;
    
    chdir "$Analysispath";
    
    #renames the table file so there is a full list of bams that can be backed up
    `cp $RNA_table $RNA_overlap` unless (-e "$RNA_overlap");
    `cp $WGS_table $WGS_overlap` unless (-e "$WGS_overlap");
    `cp $Geno_table $Geno_overlap` unless (-e "$Geno_overlap");
  
   open (GENI,"$Geno_overlap");
   open (GN,">$cancer_type\_Normal_Gen.txt");
   open (GT,">$cancer_type\_Tumor_Gen.txt");
   
    #Format of Genotypes file
    
    #UUID: 1e2b5cf3-a939-4dea-9b77-68faaf871c5a
    #TCGA ID with Sample ID: TCGA-97-7546-01A
    #Sample ID: 1
    #Sample Type: Primary Tumor
    #TCGA ID: TCGA-97-7546
    
    while (my $r = <GENI>)
    {
        chomp($r);
        my @TN = split("\t",$r);
        if ($TN[2] < 10)
        {
            print GT $r,"\n";
        }
        elsif ($TN[2] > 9)
        {
            print GN $r,"\n";
        }
    }
    close (GENI);
    close (GN);
    close (GT);
    
    #Intersect the normal and tumor Genotypes to get matching tumor normal ids
    Intersect_Files("$cancer_type\_Normal_Gen.txt,$cancer_type\_Tumor_Gen.txt","4,4","$cancer_type.Genotype_Int.txt");
    
    `rm $cancer_type\_Normal_Gen.txt $cancer_type\_Tumor_Gen.txt`;
    
    #parse files to remove non matching tumor/normal
    open (GIN,"$cancer_type.Genotype_Int.txt");
    open (GOUT,">$cancer_type\_Genotype.info.txt");
 
    while (my $r = <GIN>)
    {
        chomp($r);
        my @Genos = split("\t",$r);
        
        if ($Genos[1] =~ /normal/ig && $Genos[1] =~ /tumor/ig)
        {
            print GOUT $r,"\n";
        }
    }
    close (GIN);
    close (GOUT);
    
    #Intersect parsed file with RNA-Seq and WGS
    Intersect_Files("$RNA_table,$WGS_table,$cancer_type\_Genotype.info.txt","1,1,2","$cancer_type\_Intersect.txt");
    
    `rm $cancer_type.Genotype_Int.txt $cancer_type\_Genotype.info.txt`;
    
    #parse intersected file to remove IDs that do not intersect with RNA-Seq, WGS and Genotypes
    open (INT,"$cancer_type\_Intersect.txt");
    open (INTOUT,">$Intersect");
    
    while (my $r = <INT>)
    {
        chomp($r);
        my @RWG = split("\t",$r);
        
        if ($RWG[1] =~ /RNA/ig && $RWG[1] =~ /WGS/ig && $RWG[1] =~ /Genotype/ig)
        {
            print INTOUT $r,"\n";
        }
    }
    close (INT);
    close (INTOUT);
    
    `rm $cancer_type\_Intersect.txt`;
    
    #determine whether $Exp_Strategy is RNA-Seq of WGS
    if ($Exp_Strategy eq "RNA-Seq")
    {
	vlookup("$RNA_overlap",2,"$Intersect",3,3,"y","$Intersect\_NaN");
	`grep -v -i NaN $Intersect\_NaN > $Intersect\_pull`;
	pull_column("$Intersect\_pull","1,2,3,4,5,6,7","$RNA_overlap");
    }
    else
    {
	vlookup("$WGS_overlap",2,"$Intersect",3,3,"y","$Intersect\_NaN");
	`grep -v -i NaN $Intersect\_NaN > $Intersect\_pull`;
	pull_column("$Intersect\_pull","1,2,3,4,5,6,7,8","$WGS_overlap");
    }
    
    `rm $Intersect\_NaN $Intersect\_pull`;
}

sub Intersect_Files
{
    my $files = shift;
    my $self = $files and $files = shift if ref $files;
    my @filename = split(",",$files);
    my $IntersectColsInEachFile = shift;
    my $intersect_file = shift;
    my @cols = split(",",$IntersectColsInEachFile);
    
    my @AllItems = ();
    my %Letters_hash = ();
    my %count = ();
    my %inv_count = ();
    
    %Letters_hash=map {
	$filename[$_] => Letters($_)
    } 0..$#filename;
    my %filename_hash = reverse %Letters_hash;
    
    map
    {
	#open_and_transform;
	my $f = $filename[$_];
	my $col = $cols[$_];
	open (my $fileHandle,"$f");
	my @array;
	my $n = @array = <$fileHandle>;
	foreach my $line (@array)
	{
	    chomp $line;
	    next if $line =~ /^\s*$/;
	    my @new_array = split("\t",$line);
	    push @AllItems,$Letters_hash{$f}."_".$new_array[$col];              
	}
	close ($fileHandle); 
    }0..$#filename;
    
    foreach my $key (@AllItems)
    {
	$count{substr($key,2,)} .= substr($key,0,1);
    }
    my @tmp;
    foreach my $hash_key (keys %count)
    {
	my $group = $count{$hash_key};
	if ($inv_count{$group})
	{
	    @tmp = ($inv_count{$group});
	}
	else
	{
	    @tmp = ();
	}
	push (@tmp,$hash_key);
	$inv_count{$count{$hash_key}} = join("\n",@tmp);	
    }
    
    my %Seen = ();
    open (INF,">$intersect_file");
    foreach my $inv_key (sort keys %inv_count)
    {
	#this code is for removing duplicating targets in the @targets;
	my @targets = split("\n",$inv_count{$inv_key});
	my @unique_items = grep {!$Seen{$_}++} @targets;
	my $n = @unique_items;
	my $unique_items = join("\n",@unique_items);
	my @ls = split("",$inv_key);
	my %uniq;
	@ls = grep{!$uniq{$_}++}@ls;
	
	#Long format data for other software;
	foreach my $e (@unique_items)
	{	
	    print INF $inv_key,"\t";
	    print INF join("|",@filename_hash{sort @ls}),"\t",$e,"\n";
	}
    }    
    close (INF);
}

################Intersect_Files sub###############
sub Letters
{
    my @letters = ('a'..'z','A'..'Z');
    my $letters_num = shift;
    $letters_num = $letters_num;
    return ($letters[$letters_num]);
}

###############End of Intersect_Files sub##########

sub archive_files
{
    my $tar = Archive::Tar->new;
    my $self = shift;
    my (@files_to_archive) = @_;
    my $compress_archive = pop(@files_to_archive); #last element in the array is the name of the archive

    for (my $i = 0;$i < scalar(@files_to_archive);$i++)
    {
	my @file_list = `cat $files_to_archive[$i]`;
	foreach my $files (@file_list)
	{
	    chomp($files);
	    if ($files =~ /\t/)
	    {
		my ($file_path,$filename) = split("\t",$files); #file contains the path to files and the filenames themselves
		chdir($file_path);
		$tar->add_files($filename);
	    }
	    else
	    {
		$tar->add_files($files);
	    }
	}
    }
    $tar->write("$compress_archive");
}

1;
