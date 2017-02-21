#!/usr/bin/perl -w

package TCGA_Lib::Imputation_Plink;

use FindBin qw($Bin);
use strict;
use warnings;
use File::Copy qw(copy);
use lib "$Bin/";
use Parsing_Routines;
use FileHandle;
use MCE::Map;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

my $parsing = TCGA_Lib::Parsing_Routines->new;

###############################################################################
# Imputation_Plink.pm version 1
# Written by Tomas Babak, Zhongshan Cheng, Michael Vermeulen and Micheal Rollins-Green
#  in the Department of Biology at Queen's University.
#
# Imputation_Plink
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

sub CheckPlink{
    my $P = `which plink`;
    if($P =~ /^no plink/i) {
        print STDERR "You need to install plink and add its path in bashrc or bash_profile\n";
        exit;
    }
}

sub ParseAnnoAffxSNP6
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(PAR,$infile) or die("Can't open file $infile for input: $!");
    open(ANNO,">$outfile") or die("Can't open file $outfile for output: $!");
    while(my $r = <PAR> )
    {
        if($r =~ /^\#/)
        {
            next;
        }
        chomp($r);
        $r =~ s/\"//g;
        my @a = split(",",$r);
        if($a[4] eq '-')
        {
            $a[4] = '+';
            $a[8] =~ tr/ACGT/TGCA/;
            $a[9] =~ tr/ACGT/TGCA/;    
        }
        print ANNO join("\t",@a[0,2..4,8,9,1]),"\n";
    }
    close(PAR);
    close(ANNO);
}

sub PrepareAffxAlleleAccording2ImputationAlleles
{
    my $path = shift;
    my $self = $path and $path = shift if ref $path;
    my $database_path = shift;
    
    my $hap_tag = "nomono";
    #nomono=>ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz
    #nosing=>ALL.integrated_phase1_SHAPEIT_16-06-14.nosing.tgz 
    my $AffxInput = "$path/snp6.anno.txt";
    
    my $OneKG_legend_Path = "$database_path/ALL.integrated_phase1_SHAPEIT_16-06-14.$hap_tag";
    
    #legend file contents:
    #column 1 (id) - variant ID
    #column 2 (position) - base pair position
    #column 3 (a0) - allele labeled '0' in .hap file. This is the reference allele.
    #column 4 (a1) - allele labeled '1' in .hap file
    #column 5 (type) - SNP/INDEL/SV denotes type of biallelic variant
    #column 6 (source) - LOWCOV/EXOME denotes the type of sequencing used to generate the data at the site
    #column 7-10 (afr.aaf amr.aaf asn.aaf eur.aaf) alternate (non-reference) allele frequencies in 4 ancestral groups
    #column 11-14 (afr.maf amr.maf asn.maf eur.maf) minor allele frequencies (MAF) in 4 ancestral groups
    
    my $OUT = "$path/snp6.cd.txt";
    my %AFFX_SNPs;
    open AFFX, "$AffxInput" or die "Couldn\'t open the AFFX file: $!";
    while(<AFFX>)
    {
        my $line = $_;
        chomp($line);
        my @tp = split("\t",$line);
        if(uc($tp[1]) eq "X")
        {
            $tp[1] = 23;
        }
        
        #Exclude none number chromosomes;
        #Create hash and use chromosome and chr_pos as keys
        #Later, for indels, it is necessary to add another tag,
        #such as chr_pos_indel as the second keys of the hash
        
        #N.B.
        #The program IMPUTE2 identifies SNPs internally by their base pair positions
        #(which means that duplicate SNPs at a single position can cause problems).
        
        if($tp[1] =~ /\d+/)
        {
            if(defined($AFFX_SNPs{$tp[1]}{$tp[1].":".$tp[2].":"."SNP"}))
            {
                #For duplicate chr:pos in SNPARRAY;
                @{$AFFX_SNPs{$tp[1]}{$tp[1].":".$tp[2].":"."SNP_dup"}} = @tp[2,4,5,0,6];
            }
            else
            {	
               @{$AFFX_SNPs{$tp[1]}{$tp[1].":".$tp[2].":"."SNP"}} = @tp[2,4,5,0,6];
            }
        }
    }
    close AFFX;
    
    my %IntersectSNPs;
    
    open(OUT,">$OUT") or die "Could not output data into the file: $!\n";
    foreach my $k (sort keys %AFFX_SNPs)
    {
        if($k > 0 and $k < 24)
        {
            if($k == 23)
            {
                #ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend.gz
                open ONEKG,"zcat $OneKG_legend_Path/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend.gz |"
                or die "Couldn\'t open the OneKG reference maps: $!";	
            }
            else
            {
                open ONEKG,"zcat $OneKG_legend_Path/ALL.chr$k.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.$hap_tag.legend.gz |"
                or die "Couldn\'t open the OneKG reference maps: $!";	
            }
           
            while(<ONEKG>)
            {
                chomp($_);
                if($. > 1)
                {
                    my $l = $_;
                    my @temp = split("\t| ",$l);
                    @temp = @temp[0..4];
                    my $chrpos = $k.":".$temp[1].":".$temp[4];#chr:pos:SNP
                    if(exists($AFFX_SNPs{$k}{$chrpos}))
                    {
                        #check existance of snps=>chr:pos:SNP
                        $IntersectSNPs{$chrpos} = 1;#Recording available snps in OneKG ref db;
                        my @AffxAB = sort @{$AFFX_SNPs{$k}{$chrpos}}[1,2];
                        my @OneKGAB = sort @temp[2,3];
                        my @TranAffxAB = map
                        {
                            my $e = $_;
                            $e =~ tr/ACGT/TGCA/;
                            $e;
                        }@AffxAB;
                        @TranAffxAB = sort @TranAffxAB;
                        
                        #Because Impute2 OneKG ref|alt alleles are all in positive strand.
                        #Thus we can simply use the following conditions to tranform alleles;
                        #Otherwise, it is prone to error when when the strand info is absent!
                        
                        if($AffxAB[0] eq $OneKGAB[0] and $AffxAB[1] eq $OneKGAB[1])
                        {
                            print OUT "chr".$k, "\t", $AFFX_SNPs{$k}{$chrpos}->[0]-1, "\t",
                            $AFFX_SNPs{$k}{$chrpos}->[0], "\t", $AFFX_SNPs{$k}{$chrpos}->[3], "\t",
                            $temp[2]."|".$temp[3], "\t",$AFFX_SNPs{$k}{$chrpos}->[1], "\t",
                            $AFFX_SNPs{$k}{$chrpos}->[2],"\t",$temp[0],"\n";
                        }
                        elsif($TranAffxAB[0] eq $OneKGAB[0] and $TranAffxAB[1] eq $OneKGAB[1])
                        {
                            print OUT "chr".$k, "\t", $AFFX_SNPs{$k}{$chrpos}->[0]-1, "\t",
                            $AFFX_SNPs{$k}{$chrpos}->[0], "\t", $AFFX_SNPs{$k}{$chrpos}->[3],
                            "\t", $temp[2]."|".$temp[3], "\t",$TranAffxAB[0], "\t",
                            $TranAffxAB[1],"\t",$temp[0],"\n";
                        }
                        else
                        {
                            print STDERR "Conflict pattern for the SNP $temp[0] or  $AFFX_SNPs{$k}{$chrpos}->[3] ($chrpos) is observed as: \n";
                            print STDERR "Query Alleles: ", join("|",@AffxAB),"\n";
                            print STDERR "OneKG Alleles: ", join("|",@OneKGAB), "\n";
                            #After checking, we want that these SNPs having more than two observed alleles;
                            #In the ref, it contains two of these alleles, and only one is observed in two
                            #alleles of query SNP
                            print STDERR "\nWe decide to not keep these SNPs in their original format!\n\n";
                        }
                        
                        #For duplicate chr:pos SNPs;
                        if(exists($AFFX_SNPs{$k}{$chrpos."_dup"}))
                        {
                            $chrpos .= "_dup";
                            $IntersectSNPs{$chrpos} = 1;#Recording available snps in OneKG ref db;
                        }	    
                    }
                    else
                    {
                        next;
                    }
                }
            }
            close ONEKG;
        }
    }
    #Keep left SNPs in their original format!
    #Although these SNPs were not found in OneKGImputation Ref legend file, it will be useful!;
    foreach my $chr (keys %AFFX_SNPs)
    {
        foreach my $SNP (keys %{$AFFX_SNPs{$chr}})
        {
            unless (defined($IntersectSNPs{$SNP}))
            {
                print OUT "chr".$chr, "\t", $AFFX_SNPs{$chr}{$SNP}->[0]-1, "\t",
                $AFFX_SNPs{$chr}{$SNP}->[0], "\t", $AFFX_SNPs{$chr}{$SNP}->[3], "\t",
                $AFFX_SNPs{$chr}{$SNP}->[1]."|".$AFFX_SNPs{$chr}{$SNP}->[2],
                "\t",$AFFX_SNPs{$chr}{$SNP}->[1], "\t",$AFFX_SNPs{$chr}{$SNP}->[2],"\t",$AFFX_SNPs{$chr}{$SNP}->[4],"\n";
            }
        }
    }
    close OUT;
}

sub Prepare_Genos
{
    my $birdseed = shift;
    my $self = $birdseed and $birdseed = shift if ref $birdseed;
    my $BirdseedPath = shift;
    my $AffxAnno = shift;
    my $bases = shift;
    my $type = shift;
    
    my @a = split("\\.",$birdseed);
    
    unless (defined($type))
    {
        #Default type is for tumor/normal;
        print "You didn\'t choose which type of TCGA samples (tumor or normal) to be analyzed.\n";
        print "The default samples used for analysis will be normal samples only!\n";
        $type=2;#0 for normal, 1 for tumor, and 2 for tumor and normal;
    }
    
    my $pid = $a[1];
    my $tag = [split("-",$pid)]->[-1];
    $tag =~ s/\D+//g;
    my $tmp_f1;
    my $tmp_f2;
    my $tmp_num1 = $parsing->rand_num;
    my $tmp_num2 = $parsing->rand_num;
    my $tmp_num3 = $parsing->rand_num;
    my $tmp_num4 = $parsing->rand_num;
    $tmp_f1 = "tmp$tmp_num1".$tmp_num2;
    $tmp_f2 = "tmp$tmp_num3".$tmp_num4;
    my $AffxAnno1 = "$AffxAnno.TMP.$tmp_num1.$tmp_num3";
    `cp $AffxAnno $AffxAnno1`;
    if($tag =~ /[0-9]/)
    {
        if($type == 2)
        {

            print STDERR "working on tumor/normal sample: $pid\n";
            $parsing->strip_head("$BirdseedPath/$birdseed","$bases/$pid.t",2);
            $parsing->vlookup("$bases/$pid.t",1,"$AffxAnno1",4,"6,7,1","y","$tmp_f1");
            mk_genos("$tmp_f1","$tmp_f2");
            $parsing->pull_column($tmp_f2,"1,7,6","$bases/$pid");
        }
        elsif($type == 0)
        {
            if($tag > 9)
            {
                print STDERR "working on normal sample only: $pid\n";
                $parsing->strip_head("$BirdseedPath/$birdseed","$bases/$pid.t",2);
                $parsing->vlookup("$bases/$pid.t",1,"$AffxAnno1",4,"6,7,1","y","$tmp_f1");
                mk_genos("$tmp_f1","$tmp_f2");
                $parsing->pull_column($tmp_f2,"1,7,6","$bases/$pid");
            }
        }
        elsif($type == 1)
        {
            if($tag < 10)
            { 
                print STDERR "working on tumor sample only: $pid\n";
                $parsing->strip_head("$BirdseedPath/$birdseed","$bases/$pid.t",2);
                $parsing->vlookup("$bases/$pid.t",1,"$AffxAnno1",4,"6,7,1","y","$tmp_f1");
                mk_genos("$tmp_f1","$tmp_f2");
                $parsing->pull_column($tmp_f2,"1,7,6","$bases/$pid");
            }
        }
        else
        {
            print "Don\'t know your input sample type $type for analysis!\n" and die;
        }
    }
    if(defined $tmp_f1 || defined $tmp_f2)
    {
        `rm -f $tmp_f1 $tmp_f2 $AffxAnno1`;
    }
}

sub mk_genos
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $outfile = shift;
    
    open(MK,$infile) or die("Can't open file for input: $!\n");
    open(GEN,">$outfile")or die("Can't open file for output: $!\n");
    while(my $r = <MK> )
    {
        chomp($r);
        my @a = split("\t",$r);
        if($a[1] == 0)
        {
            print GEN $r, "\t", $a[3].",".$a[3], "\t";
        }
        elsif($a[1] == 1)
        {
            print GEN $r, "\t", $a[3].",".$a[4], "\t";
        }
        elsif($a[1] == 2)
        {
            print GEN $r, "\t", $a[4].",".$a[4], "\t";
        }
        else
        {
            print GEN $r, "\t", "0,0", "\t";
        }
        print GEN $a[5], "\n";
    }
    close(MK);
    close(GEN);
}

sub mk_map
{
    my $filename = shift;
    my $self = $filename and $filename = shift if ref $filename;
    my $RNA_map_path = shift;
    
    open(MKM,$filename) or die("Can't open $filename for input: $!");
    my $old = <MKM>;
    chomp($old);
    my @o = split("\t",$old);
    $o[0] =~ s/chr//;
    open(OUT,">$RNA_map_path/$o[0].map") or die "no out\n";
    
    while(my $r = <MKM>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        $a[0] =~ s/chr//;
        print OUT $o[0], "\t", $o[2], "\t0\t", $o[1], "\n";
        if($o[0] ne $a[0])
        {
            close OUT; open(OUT,">$RNA_map_path/$a[0].map") or die "no out\n";
        }
        @o = @a;
    }
    print OUT $o[0], "\t", $o[2], "\t0\t", $o[1], "\n";
    close(MKM);
    close(OUT);
}

sub mk_ped_list_4_plink
{
    my $chr_file = shift;
    my $self = $chr_file and $chr_file = shift if ref $chr_file;
    my $bases_dir = shift;
    my $pedlist = shift;
    my $outdir = shift;
    
    my $chr_fh=FileHandle->new("$chr_file");
    my @chrs = <$chr_fh>;
    $chr_fh->close;
    my $outfile = FileHandle->new(">$pedlist");
    my @samples = $parsing->get_only_files_in_dir("$bases_dir");
    @samples = grep{!/^\.$/}@samples;
    foreach my $r(@chrs) 
    {
        chomp($r);
        my $chr = $r;
        $chr =~ s/\.map//;
        for(my $i = 0;$i < scalar(@samples);$i++)
        {
            print STDERR "working on $chr for sample $samples[$i].\n";
            print $outfile "$chr\t$samples[$i]\t$i\t$outdir\n";
        }
    }
    $outfile->close;
}

sub Make_Plink_Bed
{
    my $chr = shift;#$chr;
    my $self = $chr and $chr = shift if ref $chr;
    my $sample = shift;#$samples[$i];
    my $sampleid = shift;#$i;
    my $path = shift;#$dir_path;
    my $bases = shift;
    my $map_dir = shift;
    
    #Copy map file to avoid of reading errors when using qsub;
    copy "$path/$map_dir/$chr.map", "$path/$map_dir/$sample.chr$chr.map";
    
    $parsing->vlookup("$path/$map_dir/$sample.chr$chr.map",2,"$bases/$sample",1,2,"y","$path/$map_dir/$sample.chr$chr.pro");
    
    ped_it_4_plink("$sample",$sampleid,"$path/$map_dir/$sample.chr$chr.pro","$path/$map_dir/$sample.chr$chr.ped");
   
    #Transform ped into BED of plink;
    print STDERR "Now making plink bed for $sample.chr$chr\n";
    `plink --map $path/$map_dir/$sample.chr$chr.map --ped $path/$map_dir/$sample.chr$chr.ped --make-bed --out $path/$map_dir/$sample.chr$chr`;
    
    unlink("$path/$map_dir/$sample.chr$chr.map",
           "$path/$map_dir/$sample.chr$chr.ped",
           "$path/$map_dir/$sample.chr$chr.pro",
           "$path/$map_dir/$sample.chr$chr.txt");
}

sub ped_it_4_plink
{
    my $sample = shift;
    my $fam = shift;
    my $infile = shift;
    my $outfile = shift;
    
    my $fh=FileHandle->new("<$infile");
    
    $fam = $fam + 1;
    
    if(-s $infile)
    {
        open(P4P,">$outfile") or die("Can't open file for output: $!");
        print P4P "FAM".$fam, "\t", $sample, "\t", 0, "\t", 0, "\t", 2, "\t", 0;
       
        while(my $r = <$fh>)
        {
            chomp($r);
            my @a = split("\t",$r);
            my @w = split(",",$a[4]);
            
            unless($w[0] =~ /^[ACGT]$/)
            {
                $w[0] = 0;
            }
            unless($w[1] =~ /^[ACGT]$/)
            {
                $w[0] = 0;
            }
            
            print P4P "\t".join("\t",@w);
        }
        print "\n";
        $fh->close;
        close(P4P);
    }
}

sub submit_shapeit
{
    my $Ref_Hap_Path = shift;
    my $self = $Ref_Hap_Path and $Ref_Hap_Path = shift if ref $Ref_Hap_Path;
    my $rna_path = shift;
    my $phased = shift;
    my $map_dir = shift;
    my $ped_dir = shift;
    my $logs = shift;
    my @cmds;
    
    $Ref_Hap_Path =~ s/\/$//;
    
    for(my $i = 1;$i < 23;$i++)
    {
        print "Working on $i for shapeit\n";
        push @cmds,"shapeit --input-ped $rna_path/$ped_dir/$i.ped $rna_path/$map_dir/$i.map --input-map $Ref_Hap_Path/genetic_map_chr$i\_combined_b37.txt --output-max $phased/$i.phased.haps $phased/$i.phased.sample --output-log $rna_path/$logs/$i.shapeit_log";
    }
    my $i = 'X';
    print "Working on $i for shapeit\n";
    push @cmds,"shapeit --input-ped $rna_path/$ped_dir/23.ped $rna_path/$map_dir/23.map --input-map $Ref_Hap_Path/genetic_map_chr$i\_nonPAR_combined_b37.txt --output-max $phased/$i.phased.haps $phased/$i.phased.sample --chrX --output-log $rna_path/$logs/$i.shapeit_log";
    
    return @cmds;
}

sub submit_all
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $OneKG_Ref_Path = shift;
    my $Phased_hap = shift;
    my $Phased_out_path = shift;
    my @cmds;
    
    open(SA,$infile) or die("Can't open file for input");
    
    while(my $r = <SA>) 
    {
        chomp($r);
        my @a = split("\t",$r);
        $a[0] =~ s/chr//;
        push @cmds,submit_impute_chr($a[0],$a[1],5e6,$OneKG_Ref_Path,$Phased_hap,$Phased_out_path);
    }
    return @cmds;    
}

##################submit all sub#######################
sub submit_impute_chr
{
    my $chr = shift;#chr_number
    my $chr_len = shift;#chr_number
    my $ww = shift;#window_size
    my $OneKG_Ref_Path = shift;
    my $Phased_hap = shift;
    my $Phased_out_path = shift;
    
    $chr =~ s/chr//i;
    
    my $impute2_m = "$OneKG_Ref_Path/genetic_map_chr$chr\_combined_b37.txt";
    my $impute2_h = "$OneKG_Ref_Path/ALL.chr$chr.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz";
    my $impute2_l = "$OneKG_Ref_Path/ALL.chr$chr.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz";
    my $know_haps_g = "$Phased_hap/$chr.phased.haps";
    my $phased_chunk_out = "$Phased_out_path/chr$chr.chunk";
    
    if($chr =~ /X/i or $chr == 23)
    {
        $impute2_m = "$OneKG_Ref_Path/genetic_map_chrX_nonPAR_combined_b37.txt";
        $impute2_h = "$OneKG_Ref_Path/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap.gz";
        $impute2_l = "$OneKG_Ref_Path/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend.gz";
        $know_haps_g = "$Phased_hap/X.phased.haps";
        $phased_chunk_out = "$Phased_out_path/chr23.chunk";
    }
    my $chunk = 1;
    my $i = 1;
    my @cmds;
        
    for( $i = 1;$i < ($chr_len - $ww);$i += $ww)
    {
        my $ss = $i;
        my $ee = $i + $ww - 1;
        push @cmds,"impute2 -use_prephased_g -m $impute2_m -h $impute2_h -l $impute2_l -known_haps_g $know_haps_g -int $ss $ee -o $phased_chunk_out.$chunk -phase";
        
        $chunk++;
    }
    my $ss = $i;
    my $ee = $chr_len;
    push @cmds,"impute2 -use_prephased_g -m $impute2_m -h $impute2_h -l $impute2_l -known_haps_g $know_haps_g -int $ss $ee -o $phased_chunk_out.$chunk -phase";
    print "We are going to do imputation for chr$chr!\n";
    return @cmds;
}
##################end of submit all sub#######################

sub impute_2_out
{
    my $path = shift;
    my $self = $path and $path = shift if ref $path;
    my $RNA_Path = shift;
    my $infile = shift;
    my $phased = shift;
    my $imputation = shift;
    
    open(IM2O,$infile) or die("Can't open file for output: $!");
    my @cmds;

    while(<IM2O>)
    {
        chomp;
        $_ =~ s/_haps//;
        my $id = $_;
        my ($chr) = $id =~ /chr(\d+)\./;
        system("cp $phased/21.phased.sample $imputation/$id.sample");
        push @cmds,"plink --oxford-single-chr $chr --hard-call-threshold 0.05 --gen $path/$id --sample $imputation/$id.sample --make-bed --out $imputation/$id";
    }
    return @cmds;
}

sub Assign_ChrPos_MissingId_SNPs
{
    my $bimpath = shift;
    my $self = $bimpath and $bimpath = shift if ref $bimpath;
    
    my %UniqSNPs;#Used to change depulicated SNP names!
    my $fh = FileHandle->new("< $bimpath") or die "Can not find the bim file: $!\n";
    my @SNPs = <$fh>;
    $fh->close;
    my $outfh = FileHandle->new(">$bimpath");
    map
    {
        my $e = $SNPs[$_];
        chomp($e);
        $e =~ s/^(\S+)\t\.\t0\t(\d+)(.*)$/$1\tchr$1:$2\t0\t$2$3/ if $e=~/^\S+\t\.\t0/;    
        my @as = split("\t| ",$e);
        if(defined($UniqSNPs{$as[1]}))
        {
            $as[1].="_".$as[3]."_".$as[4];
            print $outfh join("\t",@as),"\n";
        }
        else
        {
            $UniqSNPs{$as[1]} = 1;
            print $outfh $e,"\n";
        }
    }(0..$#SNPs);
    $outfh->close;
}

sub Extract_Ind_Het_Genos
{
    my $chr = shift;#should be "chr\d+";
    my $self = $chr and $chr = shift if ref $chr;
    my $BED = shift;
    my $phased_sample_info = shift;
    my $cds_plink_dir = shift;
    
    my $fam_fh=FileHandle->new("$BED.fam");
    chomp(my @ids=<$fam_fh>);
    $fam_fh->close;
    my @TCGAs;
    my @cmds = map
    {
        my $l = $_;
        my @es = split("\t| ",$l);
        push @TCGAs,$es[1];
        my $keep_fh = FileHandle->new(">$es[1].$chr.txt");
        print $keep_fh join("\t",@es[0..1]);
        $keep_fh->close;
         my $cmd = "plink --bfile $BED --keep $cds_plink_dir/$es[1].$chr.txt --make-bed --out $cds_plink_dir/temp.$chr.$es[1]; plink --bfile $cds_plink_dir/temp.$chr.$es[1] --max-mac 1 --mac 1 --write-snplist --out $cds_plink_dir/$es[1].$chr; rm $cds_plink_dir/temp.$chr.$es[1].*"; 
    }@ids;
    return @cmds;
}

sub Append_snplist_chrs
{
    my $ID = shift;
    my $self = $ID and $ID = shift if ref $ID;
    my $outdir = shift;
    
    my @files = glob "$ID.chr*.snplist";#In linux, if included * in the perl command line, perl script will not work!
    $outdir =~ s/\/$//;
    
    open OUT,">$outdir/$ID.snplist" or die "Couldn\'t write data into the output file $outdir/$ID.snplist: $!";
    
    my @filenames = sort
    {
        $a cmp $b
    }@files;#Important for downstream analysis, keep the order!;                 
    
    foreach my $file(@filenames)
    {
        open IN, "$file" or die "Cannot open $file for read: $!";
        while(my $l = <IN>)
        {
            chomp($l);
            print OUT $l,"\n" unless $l=~/^\s+$/;#remove empty line!;
        }
        close IN;
    }
    close OUT;
}

sub merge_tcga_snplist
{
    my $RNA_cds_plink_path = shift;
    my $self = $RNA_cds_plink_path and $RNA_cds_plink_path = shift if ref $RNA_cds_plink_path;

    for(my $chr = 1;$chr < 24;$chr++)
    {
        Get_Uniq_SNPS_From_SNPlists( "$RNA_cds_plink_path","chr$chr.snplist","Chr$chr\_ALL_SNPLIST.txt");
    }
}

sub Get_Uniq_SNPS_From_SNPlists
{
    my $dir = shift;#Is used to input and output files
    my $self = $dir and $dir = shift if ref $dir;
    my $file_tag = shift;
    my $out_filename = shift;
    
    $dir =~ s/\/$//;
   
    opendir(DIR,$dir) or die "can not open the $dir directory: $!\n";
    my @files = readdir(DIR);
    close(DIR);
       @files = grep{/$file_tag/;}@files;
    my %Uniq_line;
    open(OUT,">$dir/$out_filename") or die "Can not write unique lines into the file: $!\n";
    foreach my $f (@files)
    {
        my $FH=FileHandle->new("$dir/$f");
        while(my $l = <$FH>)
        {
            chomp($l);
            if(($l !~ /^\s*$/) and (!defined($Uniq_line{$l})))
            {
                $Uniq_line{$l} = 1;
                print OUT $l,"\n";
            }
        }
        $FH->close;
    }
    close OUT;
}

sub Snplist_Query_Haps_for_Beds
{
    my $SNPList = shift;
    my $self = $SNPList and $SNPList = shift if ref $SNPList;
    my $phased_sample_info = shift;
    my $hap_dir = shift;
    my $out_bed = shift;
    my $chr = shift;
    
    $hap_dir =~ s/\/$//;

    open(SNPList,"$SNPList")
    or die "Could not find the SNPList $SNPList: $!\n";
    my %SNPList_info;
    #Load SNPList into hash;
    while(my $l=<SNPList>)
    {
        chomp($l);
        my @ln = split("\t| ",$l);
        if($ln[1] == $chr)
        {
            my $SNP = $ln[0];
            $SNPList_info{$SNP} = 1;
        }
    }
    close SNPList;
    # load in samples
    # which contains all the sample ids;
    open(SS,$phased_sample_info) or die "no samples file\n";
    my @ss;
    my $r = <SS>;
    $r = <SS>;
    #remove first two lines;
    my $cc = 0;
    while($r = <SS>)
    {
        chomp($r);
        my @a = split("\t| ",$r);
        $ss[$cc] = $a[1];
        $cc++;
    }
    close(SS);
    
    opendir(my $HAP,$hap_dir) or die "can not open the dir $hap_dir: $!\n";
    my @HapFiles = readdir($HAP);
       @HapFiles = grep{/^chr$chr.*_haps$/;}@HapFiles;
    close($HAP);
    my $Out_FH = FileHandle->new(">$out_bed.chr$chr.bed") or die "Could not write info into $out_bed.bed: $!\n";
    print $Out_FH "Chr\tSt\tEnd\trsID\t",join("\t",@ss),"\n";
    # go through haps and geno_probs
    foreach my $hap(@HapFiles)
    {
        open(H,"$hap_dir/$hap") or die "no haps file: $!\n";
        my ($chrom) = $hap =~ /(chr\w+)\.chunk\.\d+_haps/i;
        while(my $hh = <H>)
        {
            chomp($hh);
            my @h = split("\ ",$hh);
            my @allel;
            $allel[0] = $h[3];
            $allel[1] = $h[4];
            my $rsid = $h[1];
            #Some SNPs may be without SNP IDs;
            if($rsid =~ /^\./)
            {
                $rsid = $chrom.":".$h[2];#chr\d+:pos ID;
            }
           
            if(defined($SNPList_info{$rsid}))
            {
                my $ppos = $h[2];
                
                # only focus on single base variants
                
                my $st = $ppos - 1;
                my $out = join("\t",$chrom,$st,$ppos,$rsid);
                print $Out_FH $out;
                for(my $ii = 5;$ii < @h;$ii += 2)
                {
                    my $a_pos = $h[$ii];
                    my $b_pos = $h[$ii+1];               
                    my $A = $allel[$a_pos];
                    my $B = $allel[$b_pos];
                    if($A ne $B)
                    {
                        print $Out_FH "\t",$A,"|",$B;
                    }
                    else
                    {
                        print $Out_FH "\tNA|NA"; 
                    }
                }
                print $Out_FH "\n";
                #To release memory;
                delete $SNPList_info{$rsid};
            }
        }
       close H;
    }
    close $Out_FH;
}

sub Split_Bed_Haps_for_individuals
{
    my $RNA_cdsplink_path = shift;
    my $self = $RNA_cdsplink_path and $RNA_cdsplink_path = shift if ref $RNA_cdsplink_path;
    
    my @cmds;
    for(my $chr = 1;$chr < 24;$chr++)
    { 
        #Sort bed according to chr and pos, which is important for downstream analysis!
        #Because the header is started with "Chr", it will be keeped at the top even after sorting!
        push @cmds, "sort -k 1,1 -k 2,2n $RNA_cdsplink_path/All_Hets_Haps.chr$chr.bed > $RNA_cdsplink_path/All_Hets_Haps.chr$chr.sorted.bed"
    }
    return @cmds;
}

sub Split_BedHaps
{
    my $Bed_Haps=shift;
    my $self = $Bed_Haps and $Bed_Haps = shift if ref $Bed_Haps;
    my $outdir = shift;
    my $chr = shift;
    my $filehandle_num = shift;
    
    $outdir =~ s/\/$//;
    
    my $BED_FH=FileHandle->new("$Bed_Haps") or die "Could not open the Bed Hap file: $!\n";
    chomp(my $header = <$BED_FH>);
    $BED_FH->close;
    
    my @Hds = split("\t| ",$header);
    my $opened_files_num = (@Hds-4)/$filehandle_num;#Open $filehandle_num files each time;
       ($opened_files_num) = $opened_files_num =~ /^(\d+)\.\d+/;
    
    my $fn;
    
    for($fn = 1;$fn <= $opened_files_num;$fn++)
    {
        my $BEDFH = FileHandle->new("$Bed_Haps") or die "Could not open the Bed Hap file: $!\n";
        <$BEDFH>;
        my @FHs;
        for(my $i = ($fn - 1) * $filehandle_num;$i < $filehandle_num + ($fn - 1) * $filehandle_num;$i++)
        {
            my $H = FileHandle->new(">$outdir/$Hds[$i+4].chr$chr.bed")
            or die "Can not open the bed $outdir/$Hds[$i+4].chr$chr.bed for exporting: $!\n";
            push @FHs,$H;
        }
       
       while(my $l = <$BEDFH>)
        {
            chomp($l);
            my @as = split("\t",$l);
            my $chr_pos = join("\t",@as[0..2]);
            for(my $ii = ($fn - 1) * $filehandle_num;$ii < $filehandle_num + ($fn - 1) * $filehandle_num;$ii++)
            {
                my $hi = $ii - ($fn - 1) * $filehandle_num;
                if(defined($FHs[$hi]))
                {
                    my $FileHandle = $FHs[$hi];
                    #Keep only SNPs;
                    if(length($as[$ii+4]) == 3)
                    {
                        #$as[$ii+4] eq "NA|NA" is also removed!
                        #Keep only het SNPs!;
                        print $FileHandle $chr_pos,"\t",$as[$ii+4],"\t$as[3]\n";
                    }
                }
            }
        }
       $BEDFH->close or die;
       for(my $fi = 0;$fi < @FHs;$fi++)
        {
            my $FH = $FHs[$fi];
            $FH->close or die "Can not close the $FH: $!\n";
        }
    }
    
    my $BEDFH1 = FileHandle->new("$Bed_Haps") or die "Could not open the Bed Hap file: $!\n";
    <$BEDFH1>;
    my @LeftFHs;
    for(my $li = ($fn - 1) * $filehandle_num;$li < @Hds - 4;$li++)
    {
        my $HH = FileHandle->new(">$outdir/$Hds[$li+4].chr$chr.bed")
        or die "Can not open the bed $outdir/$Hds[$li+4].chr$chr.bed for exporting: $!\n";
        push @LeftFHs,$HH;
    }
  
    while(my $l = <$BEDFH1>)
    {
        chomp($l);
        my @as=split("\t",$l);
        my $chr_pos = join("\t",@as[0..2]);
        for(my $x = ($fn - 1) * $filehandle_num;$x < @as - 4;$x++)
        {
            my $pi = $x - ($fn - 1) * $filehandle_num;
            if(defined($LeftFHs[$pi]))
            {
                my $FileHandle1 = $LeftFHs[$pi];
                #Keep only SNPs;
                if(length($as[$x+4]) == 3)
                {
                    #$as[$x+4] eq "NA|NA" is also removed!
                    #Keep only het SNPs!;
                    print $FileHandle1 $chr_pos,"\t",$as[$x+4],"\t$as[3]\n";
                }
            }
        }
    }
    
    for(my $lfi = 0;$lfi < @LeftFHs;$lfi++)
    {
        my $LeftFH = $LeftFHs[$lfi];
        $LeftFH->close or die "Can not close the $LeftFH: $!\n";    
    }
    $BEDFH1->close;
}

sub Merge_All_Chrs_4_Single_TCGA_Bed
{
    my $sample_file = shift;
    my $self = $sample_file and $sample_file = shift if ref $sample_file;
    my $snplist_dir = shift;
    my $out_dir = shift;
    
    my @fam_list = `cat $sample_file`;
    @fam_list = grep {/TCGA-*-*-*/}@fam_list;
  
    for (my $i = 0;$i < scalar(@fam_list);$i++)
    {
            chomp($fam_list[$i]);
            my @tcga_id = split("\t| ",$fam_list[$i]);
            my $id = $tcga_id[1];
            Append_Beds("$id","$snplist_dir","$out_dir");       
    }
}

sub Append_Beds
{
    my $ID = shift;
    my $self = $ID and $ID = shift if ref $ID;
    my $snplist_dir = shift; 
    my $outdir = shift;
    
    my @files = glob "$ID.chr*";#In linux, if included * in the perl command line, perl script will not work!
    $snplist_dir =~ s/\/$//;
    $outdir =~ s/\/$//;
    #Import snplist of TCGAID to remove SNPs with geno probability >0.95;
    my $SNPLIST=FileHandle->new("$snplist_dir/$ID.snplist") or die "can not open the $ID.snplist: $!\n";
    my %SNPs;
    while(my $snp = <$SNPLIST>)
    {
        chomp($snp);
        $SNPs{$snp} = 1;   
    }
    $SNPLIST->close;
    my $Out="$outdir/$ID.bed";
    open OUT,">$Out" or die "Couldn\'t write data into the output file $Out: $!";
    my @filenames = sort
    {
        $a cmp $b
    }@files;                                   
    foreach my $file(@filenames)
    {
        open IN, "$file" or die "Cannot open $file for read: $!";
        my $i = 0;
        while(my $l = <IN>)
        {
            chomp($l);
            my @as = split("\t| ",$l);
            if(defined($SNPs{$as[-1]}))
            {
                unless ($l =~ /^Chr\S+/)
                {
                    #remove header;
                    $l =~ s/^chr23/chrX/i;
                    print OUT $l,"\n";
                }
            } 
        }
        close IN;
    }
    close OUT;
}

sub Change_CD_ChrInfo_for_Mpileup
{
    my $cd = shift;
    my $self = $cd and $cd = shift if ref $cd;
    my $cd_dir = shift;
    my $outfile = shift;
    
    $cd_dir =~ s/\/$//;
    
    open(COUT,">$outfile") or die("Can't open file for output: $!");
    my $FH = FileHandle->new("$cd_dir/$cd") or die "no cd input: $!\n";
    while(my $l = <$FH>)
    {
        chomp($l);
        $l =~ s/^chr23/X/i;
        $l =~ s/^chr//i;
        print COUT $l,"\n"; 
    }
    $FH->close;
    close(COUT);
}

sub fetch_Chrom_Sizes
{
    my $DB = shift;
    my $self = $DB and $DB = shift if ref $DB;
    system '/bin/sh', '-c', <<'StuffyFunk', $0, $DB;
#!/bin/sh

# fetchChromSizes - script to grab chrom.sizes from UCSC via either of:
#	mysql, wget or ftp

# "$Id: fetchChromSizes,v 1.2 2009/03/26 18:40:09 hiram Exp $"

usage() {
    echo "usage: fetchChromSizes <db> > <db>.chrom.sizes"
    echo "   used to fetch chrom.sizes information from UCSC for the given <db>"
    echo "<db> - name of UCSC database, e.g.: hg18, mm9, etc ..."
    echo ""
    echo "This script expects to find one of the following commands:"
    echo "   wget, mysql, or ftp in order to fetch information from UCSC."
    echo "Route the output to the file <db>.chrom.sizes as indicated above."
    echo ""
    echo "Example:   fetchChromSizes hg18 > hg18.chrom.sizes"
    exit 255
}

DB=$1
export DB
if [ -z "${DB}" ]; then
    echo "ERROR: No database was specified"
    exit 255
fi

WGET=`type wget 2> /dev/null | sed -e "s/.* is //"`
MYSQL=`type mysql 2> /dev/null | sed -e "s/.* is //"`
FTP=`type ftp 2> /dev/null | sed -e "s/.* is //"`

if [ -z "${WGET}" -a -z "${MYSQL}" -a -z "${FTP}" ]; then
    echo "ERROR: can not find any of the commands: wget mysql or ftp"
    exit 255
fi

DONE=0
export DONE
if [ -n "${MYSQL}" ]; then
    echo "INFO: trying MySQL ${MYSQL} for database ${DB}" 1>&2
    ${MYSQL} -N --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
	"SELECT chrom,size FROM chromInfo ORDER BY size DESC;" > chr_lens_grep_chr ${DB} 
    if [ $? -ne 0 ]; then
	echo "WARNING: mysql command failed" 1>&2
    else
	DONE=1
    fi
fi

if [ -n "${WGET}" ]; then
    echo "INFO: trying WGET ${WGET} for database ${DB}" 1>&2
    tmpResult=chromInfoTemp.$$.gz
    ${WGET} \
"ftp://hgdownload.cse.ucsc.edu/goldenPath/${DB}/database/chromInfo.txt.gz" \
	-O ${tmpResult} 2> /dev/null
    if [ $? -ne 0 -o ! -s "${tmpResult}" ]; then
	echo "WARNING: wget command failed" 1>&2
	rm -f "${tmpResult}"
    else
	zcat ${tmpResult} 2> /dev/null | cut -f1,2 | sort -k2nr > chr_lens_grep_chr
	rm -f "${tmpResult}"
        DONE=1
    fi
fi

if [ -n "${FTP}" ]; then
    echo "INFO: trying FTP ${FTP} for database ${DB}" 1>&2
    tmpFtpRsp=fetchTemp.$$
    echo "user anonymous fetchChromSizes@" > ${tmpFtpRsp}
    echo "cd goldenPath/${DB}/database" >> ${tmpFtpRsp}
    echo "get chromInfo.txt.gz ${tmpResult}" >> ${tmpFtpRsp}
    echo "bye" >> ${tmpFtpRsp}
    ${FTP} -u -n -i hgdownload.cse.ucsc.edu < ${tmpFtpRsp} > /dev/null 2>&1
    if [ $? -ne 0 -o ! -s "${tmpResult}" ]; then
	echo "ERROR: ftp command failed" 1>&2
	rm -f "${tmpFtpRsp}" "${tmpResult}"
    else
	zcat ${tmpResult} | cut -f1,2 | sort -k2nr > chr_lens_grep_chr
	rm -f "${tmpFtpRsp}" "${tmpResult}"
        DONE=1
    fi
fi

if [ "${DONE}" -eq 0 ]; then
    echo "ERROR: sorry, attempt to fetch chrom.sizes has failed ?" 1>&2
    exit 255
fi

StuffyFunk
}

1;
