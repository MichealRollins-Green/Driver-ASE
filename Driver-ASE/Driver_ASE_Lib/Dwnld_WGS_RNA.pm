#!/usr/bin/perl -w
package Driver_ASE_Lib::Dwnld_WGS_RNA;

use strict;
use warnings;
use FileHandle;
use LWP::Simple;
use FindBin qw($Bin);
use File::Copy qw(copy);
use MCE::Map max_workers => 4;
use Cwd;
use Cwd 'realpath';
my $TCGA_Analysis_Dir = realpath("../");
use lib "$Bin/";
use Parsing_Routines;
use autodie;
no warnings 'once';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

my $parsing = Driver_ASE_Lib::Parsing_Routines->new;

###############################################################################
# Dwnld_WGS_RNA.pm version 1
# Written by Tomas Babak, Zhongshan Cheng, Michael Vermeulen and Micheal Rollins-Green
#  in the Department of Biology at Queen's University.
#
# Dwnld_WGS_RNA
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

#parses the gdc website manifest of the specified cancer type and prints to a results file.
sub gdc_parser
{
    my $cancer_type = shift;
    my $self = $cancer_type and $cancer_type = shift if ref $cancer_type;
    my ($exp_strat,$array_type) = @_; #Experimental Strategy, array type (Genotypes or Copy number estimate)
    
    my @cancer = split(" ",$cancer_type);
    
    foreach my $cancer_type(@cancer)
    {
        #Go to this web to change hex into ASCII;
        #https://en.wikipedia.org/wiki/ASCII;
        #The gdc html link should be: https://gdc-api.nci.nih.gov/legacy/files?;
        my $url;
       if ($exp_strat =~ /Genotyping array/)
        {
           $url = "https://gdc-api.nci.nih.gov/legacy/files?filters=%7B'op':'and','content':%5B%7B'op':'in','content':%7B'field':'cases.project.project_id','value':%5B'TCGA-$cancer_type'%5D%7D%7D,%7B'op':'in','content':%7B'field':'files.data_type','value':%5B'$array_type'%5D%7D%7D%5D%7D";
        }
        elsif ($exp_strat =~ /RNA-Seq/ || $exp_strat =~ /WGS/)
        {
            $url = "https://gdc-api.nci.nih.gov/legacy/files?filters=%7B'op':'and','content':%5B%7B'op':'in','content':%7B'field':'files.data_type','value':%5B'Aligned%20reads'%5D%7D%7D,%7B'op':'in','content':%7B'field':'files.experimental_strategy','value':%5B'$exp_strat'%5D%7D%7D,%7B'op':'in','content':%7B'field':'cases.project.project_id','value':%5B'TCGA-$cancer_type'%5D%7D%7D%5D%7D";
        }
        #my $url="https://dcc.icgc.org/api/v1/manifests?filters=%7B'donor':%7B'primarySite':%7B'is':%5B'Blood'%5D%7D%7D%7D";
        $url =~ s/https://;
        $url =~ s/'/%22/g;
        $url =~ s/"/%22/g;
        $url =~ s/,/%2C/g;
        $url =~ s/:/%3A/g;
        $url =~ s/\{/%7B/g;
        $url =~ s/\}/%7D/g;
        $url =~ s/ /%20/g;
        $url =~ s/\]/%5D/g;
        $url =~ s/\[/%5B/g;
        $url = "https:$url&size=30000&return_type=manifest";
        print "Now we are querying all $cancer_type cancer samples from TCGA project:\n\n",$url,"\n\n";
        if (defined $array_type)
        {
            `curl \'$url\' > \'$cancer_type.$array_type.result.txt\'`;
        }
        else
        {
            `curl \'$url\' > \'$cancer_type.result.txt\'`;
        }
    }
}

#Gets the metadata data of the file and places the UUID in a Payload.txt file.
sub metadata_collect
{
    my $result_file = shift; #result file made from gdc_parser
    my $self = $result_file and $result_file = shift if ref $result_file;
    my $payload = shift; #parsed results outputted to payload file
    open (my $in,$result_file);
    open (MQ,">$payload");
    chomp(my @uuids = <$in>);
    close ($in);
    shift @uuids;
    @uuids=grep{$_ = "\t\t".'"'.[split("\t| ", $_)]->[0].'"' unless $_ =~ /^\s*$/;}@uuids;
    print MQ '{
        "filters":{
            "op":"in",
            "content":{
                "field":"files.file_id",
                "value":[',"\n";
    my $last=pop @uuids;
    print MQ join(",\n",@uuids),",\n",$last,"\n";
    
    print   MQ '             ]
            }
        },
        "format":"TSV",
    "fields":"file_id,cases.submitter_id,cases.case_id,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,platform",
        "size":"30000"
    }',"\n";
    close (MQ);
}

#Uses the metadata file to get data(i.e. UUID, TCGA ID) and prints the data for the UUIDs to a datatable file.
sub parse_patient_id
{
    ################################################################
    # Parses the tumour-normal tissue number out of the Patient ID # 
    ################################################################
    
    my $metadata_file = shift; #metadata file
    my $self = $metadata_file and $metadata_file = shift if ref $metadata_file;
    my $datatable_file = shift; #output data to a datatable file
    open (my $in, $metadata_file); #open input_file
    chomp(my @data = <$in>);
    close $in;
    open (PID,">$datatable_file");
    shift @data;  #remove headers
    
    foreach my $lines(@data)
    {
        my @split_line = split("\t",$lines);
        my $Patient_ID_column = $split_line[1];
        $Patient_ID_column = substr $Patient_ID_column, -3, 2;

       if ($Patient_ID_column == 01)
        {
            $Patient_ID_column = 1;
        }
       if ($Patient_ID_column == 02)
        {
            $Patient_ID_column = 2;
        }
        #.dtable(file.id | Patient ID | Sample Type | Tumour/Normal(num) | Platform)
        print PID $split_line[4], "\t", $split_line[5], "\t", $Patient_ID_column, "\t", $split_line[0], "\t", $split_line[3], "\n";
    }
    close (PID);
}

sub metadata_ids
{
    my $result_file = shift; #result file from gdc_parser
    my $self = $result_file and $result_file = shift if ref $result_file;
    my $payload = shift; #output data to a payload file
    open (my $in,$result_file);
    open (PAY,">$payload");
    chomp(my @uuids = <$in>);
    close ($in);
    shift @uuids;
    @uuids = grep{$_ = "\t\t".'"'.[split("\t| ", $_)]->[0].'"' unless $_ =~ /^\s*$/;}@uuids;
    print PAY '{
        "filters":{
            "op":"in",
            "content":{
                "field":"files.file_id",
                "value":[',"\n";
    my $last=pop @uuids;
    print PAY join(",\n",@uuids),",\n",$last,"\n";
    
    print  PAY  '             ]
            }
        },
        "format":"TSV",
    "fields":"file_id,metadata_files.file_id,files.metadata_files.file_id,metadata_files.file_name",
        "size":"30000"
    }',"\n";
    close (PAY);
}

sub parse_meta_id
{
    my $meta_table=shift;
    my $self = $meta_table and $meta_table = shift if ref $meta_table;
    my ($UUID_list,$meta_ids) = @_; #file with UUIDs and file with meta IDs

    my @ID_column;
    
    open (UUIDS, "$UUID_list");
    chomp(my @all_bam_uuids = <UUIDS>);
    close UUIDS;
    
    open (my $fh,"$meta_table");
    chomp(my @table = <$fh>);
    close $fh;
    open (PMI,">$meta_ids");
    shift @table; #remove headers
    
    foreach my $line(@table)
    {
        my @split_line = split("\t", $line);
        my @ids = grep{!/\.xml/}@split_line;
        my $bam_uuid;
        map
        {
            my $id = $_;
            foreach my $i(@all_bam_uuids)
            {
                $bam_uuid = $id if ($id eq $i);
            }
        }@ids;
        map
        {
            my $xml_id = $_;
            push @ID_column,"$xml_id\t$bam_uuid" unless( $xml_id eq $bam_uuid); 
        }@ids;
    }
    
    print PMI join("\n", @ID_column);
    print PMI "\n";
    close (PMI);
}

sub ref_parse
{
    my $meta_ID = shift; #meta IDs file
    my $self = $meta_ID and $meta_ID = shift if ref $meta_ID;
    my ($gdc_key_path,$reference,$dwld_cmd) = @_; #path to key, reference file, download command
    
    my %refs;
    my $ref_tmp = "ref_tmp";
    my $token = ImportAlldataIntoVar("$gdc_key_path/gdc.key");
    
    open (my $fh,$meta_ID);
    my @table = <$fh>;
    close ($fh);
    open (RP,">$reference");
    mkdir "$ref_tmp" unless (-d "$ref_tmp");

    my $cmd;
   mce_map
   {
        chomp($_);
        my($refID,$UUID) = ($_ =~ /(\S*)\t(\S*)/);
        if ($dwld_cmd eq "curl")
        {
            $cmd = "$dwld_cmd -# -C - --header \'X-Auth-Token: $token\' \'https://gdc-api.nci.nih.gov/legacy/data/$refID\' --output \'$ref_tmp/$refID.$UUID\' ";
        }
        else
        {
            #aria2c is the download command if curl was not specified as the command in the command line when running script
            $cmd = "$dwld_cmd -s 16 -x 16 -c --header \'X-Auth-Token: $token\' \'https://gdc-api.nci.nih.gov/legacy/data/$refID\' -o \'$ref_tmp/$refID.$UUID\' ";
        }
        `$cmd`;
    }@table;
    opendir (REF, "$ref_tmp");
    my @ref_txts = grep{!/^\./ && -f "$ref_tmp/$_"} readdir(REF);
    closedir (REF);
    
    foreach my $txt(@ref_txts)
    {
        my($RID,$UID) = split("\\.",$txt);
        $/ = undef;
        open (JSON,"$ref_tmp/$txt");
        my $info = <JSON>;
        my($Nill,$ref_info) = ($info=~/<ASSEMBLY>
           (.*?)<
            STANDARD\W+short_name=\"([^\n]+)\"
            \/>
            /ismx);#STANDARD short_name=\"(.*)\"
        $refs{$UID} = $ref_info if defined($ref_info);
        close JSON;
        $/ = "\n";
    }
    `rm -rf $ref_tmp`;
    foreach my $key(keys %refs)
    {
        print RP "$key\t$refs{$key}\n";
    }
    close (RP);
}

sub index_ids
{
    my $result_file = shift; #file from gdc_parser
    my $self = $result_file and $result_file = shift if ref $result_file;
    my $payload = shift; #output to payload file
    open (my $in,$result_file);
    chomp(my @uuids = <$in>);
    close ($in);
    open (INDEX,">$payload");
    shift @uuids;
    @uuids=grep{$_ = "\t\t".'"'.[split("\t| ", $_)]->[0].'"' unless $_ =~ /^\s*$/;}@uuids;
    print INDEX '{
        "filters":{
            "op":"in",
            "content":{
                "field":"files.file_id",
                "value":[',"\n";
    my $last = pop @uuids;
    print INDEX join(",\n",@uuids),",\n",$last,"\n";
    
    print  INDEX  '             ]
            }
        },
        "format":"TSV",
    "fields":"index_files.file_id,index_files.file_name,file_id,",
        "size":"30000"
    }',"\n";
    close (INDEX);
}

sub pull_matched_tn_GDC
{
    my $sorted_downloadtable = shift; #sorted downloadtable file
    my $self = $sorted_downloadtable and $sorted_downloadtable = shift if ref $sorted_downloadtable;
    my $final_downloadtable = shift; #downloadtable that will be used to download data
    
    my $total_tn = 0;
    open (MI,"$sorted_downloadtable");
    chomp(my @lines = <MI>);
    my %ID_ref_platfm;
    my %uniq_ID_ref_platfm = map
    {
        my $l = $_;
        my @as = split("\t",$l);
        my $key = join(":",@as[1,4,5]);
           $ID_ref_platfm{$key.":".$as[3].":".$as[2]} = $l;
           $key => 1
    }@lines;
    close (MI);    
    open (MO,">$final_downloadtable");
    my @all_keys = keys %ID_ref_platfm;
    foreach my $k(keys %uniq_ID_ref_platfm)
    {
        my @targets = grep{/$k/i;}@all_keys;
       if (@targets > 1)
        {            
            my @tmp = sort{$a->[1] cmp $b->[1]}
            map
            {
                my $t = $_;
                my @xs = split(":",$t);
                my $v = pop @xs;
                [$t,$v];
            }@targets;
            foreach my $vv(@tmp)
            {
               if ($vv->[1] < 10)
                {
                    foreach my $xx(@tmp)
                    {
                       if ($xx->[1] > 9)
                        {
                            $total_tn++;
                            print MO $ID_ref_platfm{$vv->[0]},"\n";
                            print MO $ID_ref_platfm{$xx->[0]},"\n";
                        }
                    }
                }
            }
        }
    }
    close (MO);
   if ($total_tn == 0)
    {
        print STDERR "No tumor-normal pairs in file $sorted_downloadtable!\n";
        exit;
    }
}

sub Dwld_WGSBam_and_do_mpileup
{
    my $bamlist = shift; #final download file that conatains the list of bams to download
    my $self = $bamlist and $bamlist = shift if ref $bamlist;
    my ($key_dir,$wgs_output_dir,$ref_fullpath,$mpileup_outdir,$action,$pairs,$VarScan_Path,$cancer_type,$alt_ref,$tables,$dwld_cmd,$samtools) = @_; #path to gdc key, path to wgs output directory, path to reference file GRCh37-lite.fa, path to where mpileups will be stored, action (download, all or mpileups), # of pairs to download, path to VarScan jar file, cancer type, path to alternate reference file hg19.fa, path to cancer type tables directory, download command, samtools command
    
    $mpileup_outdir =~ s/\/$//;
    
   if ($action eq "all" or $action eq "download")
    {
        if (-e ("$tables/already_done_WGS.txt"))
        {
            #Checks if there a NaN in the bamlist file and removes them as it will interfere with the code below.
            open (BAMI,"$bamlist");
            open (BAMO,">$bamlist.parse");
            
            while (my $r = <BAMI>)
            {
                chomp($r);
                
                if ($r =~ /NaN/g)
                {
                    $r =~ s/NaN//;
                }
                print BAMO $r,"\n";   
            }
            close (BAMI);
            close (BAMO);
            
            #Matches the tumor bams in the bam list
            $parsing->vlookup("$bamlist.parse",1,"$tables/already_done_WGS.txt",1,1,"y","wgs_tum.txt");
            #greps all that didn't match
            `cat wgs_tum.txt|grep NaN > wgs_NaN.txt`;
            #Gets all columns except for the NaN column
            $parsing->pull_column("wgs_NaN.txt","1,2,3,4,5,6,7","wgs_get_norm.txt");
            #Matches the normal bams
            $parsing->vlookup("wgs_get_norm.txt",1,"$tables/already_done_WGS.txt",2,2,"y","wgs_norm.txt");
            #Gets the bams that did not meet any of the above criteria
            `cat wgs_norm.txt|grep NaN > $bamlist.new`;
            `rm wgs_*.txt`;
        }
        elsif (-e("$key_dir/already_done_WGS.txt"))
        {
            #Checks if there a NaN in the bamlist file and removes them as it will interfere with the code below.
            open (BAMI,"$bamlist");
            open (BAMO,">$bamlist.parse");
            
            while (my $r = <BAMI>)
            {
                chomp($r);
                
                if ($r =~ /NaN/g)
                {
                    $r =~ s/NaN//;
                }
                print BAMO $r,"\n";   
            }
            close (BAMI);
            close (BAMO);
            
            #Matches the tumor bams in the bam list
            $parsing->vlookup("$bamlist.parse",1,"$key_dir/already_done_WGS.txt",1,1,"y","wgs_tum.txt");
            #greps all that didn't match
            `cat wgs_tum.txt|grep NaN > wgs_NaN.txt`;
            #Gets all columns except for the NaN column
            $parsing->pull_column("wgs_NaN.txt","1,2,3,4,5,6,7","wgs_get_norm.txt");
            #Matches the normal bams
            $parsing->vlookup("wgs_get_norm.txt",1,"$key_dir/already_done_WGS.txt",2,2,"y","wgs_norm.txt");
            #Gets the bams that did not meet any of the above criterea
            `cat wgs_norm.txt|grep NaN > $bamlist.new`;
            `rm wgs_*.txt`;
        }
        else
        {
            `cp $bamlist "$bamlist.new"`;
        }
        
        my $newdir = $parsing->temp_filename;
        mkdir "$key_dir/$newdir" unless (-d "$key_dir/$newdir");
        chdir "$key_dir/$newdir";
        `split -l $pairs $bamlist.new`;
        chdir "$key_dir";
        opendir (CURR, "$key_dir/$newdir");
        my @fs = grep{/^x[a-z]{2}$/} readdir(CURR);
        closedir (CURR);
        
        #Going to do downloading and do mpileups for these bams.
        
        foreach my $f(sort @fs)
        {
            print "Now start to download bams in the file $f!\n";
            #dwnld_wgs_or_rna(bamlist,key directory,bam output directory)
            dwnld_wgs_or_rna("$key_dir/$newdir/$f","$key_dir","$wgs_output_dir",$dwld_cmd,"$samtools");
           
          if ($action eq "all")
            {
                my @bam_pairs = mk_files_for_wgs("$key_dir/$newdir/$f","$wgs_output_dir");
                
                launch_wgs_mpileup_and_VarScan("$ref_fullpath","$mpileup_outdir",\@bam_pairs,$wgs_output_dir,$VarScan_Path,$cancer_type,$alt_ref,"$tables","$samtools","$key_dir");
            }
        }
    }
    elsif ($action eq "mpileups")
    {
        my @bam_pairs = mk_files_for_wgs("$bamlist","$wgs_output_dir");
        
        launch_wgs_mpileup_and_VarScan("$ref_fullpath","$mpileup_outdir",\@bam_pairs,$wgs_output_dir,$VarScan_Path,$cancer_type,$alt_ref,"$tables","$samtools","$key_dir");
    }
}

sub mk_files_for_wgs
{
    my $bam_file = shift; #file that contains list of bams to download
    my $self = $bam_file and $bam_file = shift if ref $bam_file;
    my $wgs_bam_dir = shift;
    
    open (FFW,$bam_file);
    
    $wgs_bam_dir =~ s/\/$//;
    my @rst;
    while(my $r = <FFW>)
    {
        chomp($r);
        my @a = split("\t",$r);
        my $tum = $a[0];
        
        my $n = <FFW>;
        chomp($n);
        #Going to process normal sample info;
        my @na = split("\t",$n);
        my $norm=$na[0];
        
        unless($a[2] < 10 and $na[2] > 9)
        {
            print "No normal and tumor samples for the following line: \n", $n,"\n",
            $r, "\n";
            next;
        }      
        
        if (-e("$wgs_bam_dir/$tum/$tum.bam.bai") and -e("$wgs_bam_dir/$norm/$norm.bam.bai"))
        {
            my $T_bam = "$wgs_bam_dir/$tum/$tum.bam";
            my $N_bam = "$wgs_bam_dir/$norm/$norm.bam";
            
           if (-f "$T_bam" and -f "$N_bam")
            {
                push @rst, "$N_bam\t$T_bam\t$a[1]\n";
            }
        }
    }
    close FFW;
    return @rst;
}

sub launch_wgs_mpileup_and_VarScan
{
    my $fasta = shift; #path to reference file GRCh37-lite.fa
    my $self = $fasta and $fasta = shift if ref $fasta;
    my ($mpileups_out_dir,$bam_fullpaths_ref,$wgs_output_dir,$VarScan_Path,$cancer_type,$alt_ref,$tables,$samtools,$key_dir) = @_; #path to directory where mpileups will be stored, path to bam tumor/normal files, path to directory where wgs download data is stored, path to VarScan jar file, cancer type, pat to alternate reference file hg19.fa, path to cancer type table directory, samtools command, path to directory where files were processed for downloading
    
    my $tumor;
    my $normal;
    my $done_wgs; #file that will be written to when WGS bams have finish mpileups
    foreach my $tn_bam(@$bam_fullpaths_ref)
    {
    	if (-e "$tables/already_done_WGS.txt")
    	{
	    $done_wgs = "$tables/already_done_WGS.txt";
    	}
    	else
    	{
            $done_wgs = "$key_dir/already_done_WGS.txt";
    	}
        chomp($tn_bam);
        next if $tn_bam =~ /^\s*$/;
        my @a = split("\t",$tn_bam);
        
        $a[0] =~ /([^\/]+)\/[^\/]+\.bam$/;#revised!
        my $norm = $1;
        $norm .= ".".$a[2];
        print $norm," has been submitted\n";
        
        my $reference = $fasta;
        my $normal_bam = $a[0];
        my $tumor_bam = $a[1];
        
        #Check the size of bam;
        my $Gb = 1000 * 1000;
       if ((-s "$tumor_bam") < $Gb and (-s "$normal_bam") < $Gb)
        {
            print STDERR "Your bams are less than 1Gb!\n";
            print STDERR "tumor_bam: $tumor_bam\n";
            print STDERR "normal_bam: $normal_bam\n";
            exit;#stop all mpileups and need to check GDC, which may be under maintanence;
        }
	
	if (!(-e "$tumor_bam.bai") or !(-e "$normal_bam.bai"))
        {
            print "Either $tumor_bam or $normal_bam is not indexed.\n";
            print "Trying next pair.\n";
            next;
        }
        
        #Important!
        #Only with the right fasta ref for WGS Bams,
        #Varscan can give correct results;
        
        #Check compatibility between fasta ref and ref of WGS Bams;
        my $Ref_chr_label = `head -n 1 $reference`;
        print "fasta header is: \n $Ref_chr_label.\n";
        my ($chr_info) = $Ref_chr_label =~ />((chr)*\d+)\s/msi;
        my $chr_tag;
        print "This is the parsed chr info \"$chr_info\" ","for $reference.\n";
       if ($chr_info =~ /chr/i)
        {
           $chr_tag = 0;
        }
        else
        {
           $chr_tag = 1;
        }
        
        my $tum_ref_tag = Bam_Ref_Checker($tumor_bam,"$samtools");
        print "Tumor reference tag is $tum_ref_tag.\n";
        my $norm_ref_tag = Bam_Ref_Checker($normal_bam,"$samtools");
        print "Normal reference tag is $norm_ref_tag.\n";
        
        if ($chr_tag == $tum_ref_tag and $chr_tag == $norm_ref_tag)
        {
            print "The fasta reference of $reference is compatable with $tumor_bam and $normal_bam.\n";
        }
        else
        {
            #try to use alternative fasta;
            print "The fasta reference of $reference is NOT compatable with $tumor_bam and $normal_bam!\n";
            $reference = $alt_ref;
            print "Trying fasta reference of $reference, in which the chromosomes are labeled with 'chr'.\n";
        } 
        my $NT_pileup = "$samtools mpileup -q 1 -B -f $reference $normal_bam $tumor_bam";
        my $output = "$mpileups_out_dir/$norm";
        
        #Checks for existing files and goes to the next pair if the file exists. Delete files that may be incomplete.
        if (-e "$mpileups_out_dir/$norm.*.*")
        {
            next;
        }
        
        my $cmd="bash -c \"java -jar $VarScan_Path somatic <($NT_pileup) $output.varscan --mpileup 1\"";
        print "Running: \n$cmd\n\n";
        `$cmd`;
        
        $tumor = [split("/",$tumor_bam)]->[-1];
        $normal = [split("/",$normal_bam)]->[-1];
        $tumor =~ s/.bam//;
        $normal =~ s/.bam//;
	
        print "$tumor and $normal have finished.\n";
	open (DON,">>$done_wgs");
        
        print DON $tumor,"\t",$normal,"\n";
        
	close (DON);
    }
    Delete_Files_Recursively("$wgs_output_dir","\.bam");
}

######################launch_wgs_mpileup_and_VarScan sub#######################################
sub Bam_Ref_Checker
{
    my $Bam = shift;
    my $self = $Bam and $Bam = shift if ref $Bam;
    my $samtools = shift; #samtools command
    my $tag;
    unless (-e "$Bam")
    {
        print STDERR "$Bam can not be located!\n\n$!.\n";
        exit;
    }
    else
    {
        my $ls=`$samtools view -H $Bam`;
        #print $ls,"\n";
        my @as = split("\n",$ls);
        foreach my $chr_ref(@as)
        {
            my $a = $chr_ref;
            ($tag) = $a =~ /\bSN:(\S+)\s/;
            if (defined $tag)
            {
               if ($tag =~ /^chr/i)
                {
                   print "Bam is aligned to a reference with chromosome labeled as 'Chr*'\n","The corresponding reference will be used!\n";
                   print $tag,"\n";
                }
                else
                {
                   print "Bam is aligned to a reference with chromosome labeled in number\n","The corresponding reference will be used!\n";
                   print $tag,"\n";                        
                }
                last;    
            }
        }
    }
    return $tag;
}
######################end of launch_wgs_mpileup_and_VarScan sub#######################################

sub Dwld_RNASeq_Bam_and_do_mpileup
{
    my $bamlist = shift; #final download file that conatains the list of bams to download
    my $self = $bamlist and $bamlist = shift if ref $bamlist;
    my ($key_dir,$bam_output_dir,$ref_fullpath,$mpileup_outdir,$cds_sorted_path,$action,$line_num2split,$dwld_cmd,$samtools) = @_; #path to gdc key file, path to reference file GRCh37-lite.fa, path to directory where mpileups will be stored, path to sorted bed files, action (download, all or mpileups), number of bams to download at a time, download command, samtools command
    
    $bamlist =~ s/\/$//;
    print STDERR "Bam list $bamlist does not exist!\n" and exit unless(-e "$bamlist"); 
    $key_dir =~ s/\/$//;      
    print STDERR "Directory $key_dir does not exist!\n" and exit unless(-d "$key_dir");  
    $mpileup_outdir =~ s/\/$//;
    `mkdir -p $mpileup_outdir` unless(-d "$mpileup_outdir");
    
   if ($action eq "all" or $action eq "download")
    {
        if ($action eq "all")
        {
            my $BAMs;
            #Filters through the table to find the matching bed files in cds_sorted and only includes those entries that have a matching bed in the cds_sorted directory.
            `ls $cds_sorted_path > cds_beds.txt`;
            open (CDS,"cds_beds.txt");
            open (CDO,">cds_beds_out.txt");
            while(my $r = <CDS>)
            {
                chomp($r);
                $r = substr($r,0,12);
                print CDO $r, "\n";
            }
            close (CDO);
            close (CDS);
            
            $parsing->vlookup("$bamlist",2,"cds_beds_out.txt",1,1,"y","cds_beds_lookup.txt");
            `grep -v NaN cds_beds_lookup.txt > cds_beds_new.txt`;
            my $bamlist_cds = "cds_beds_new.txt";
            open ($BAMs,"$bamlist_cds");
            
            #Get already done pileup files and remove it from the $bamlist;
            print "Going to check already done pileups in the $mpileup_outdir directory.\n";
            my @already_done = `ls $mpileup_outdir|grep '.'`;
            @already_done = grep{chomp;-s "$mpileup_outdir/$_";}@already_done;
            my %already_done = map{my($k,$v) = $_ =~ /([^\.]+)\.([^\.]+)/;$k=>$v}@already_done;
            open (my $NBAMs,">$bamlist.new");
            while(my $b = <$BAMs>)
            {
                chomp($b);
                my @es = split("\t",$b);
               if (defined $already_done{$es[0]})
                {
                    print "The TCGA ID $es[0].$es[1] has pileup result in the $mpileup_outdir and will not be submitted for pileup again\n";
                }
                else
                {
                    print $NBAMs $b,"\n" ;          
                }
            }
            close ($BAMs);
            close ($NBAMs);
            `cp $bamlist $bamlist.backup -f`;
            $bamlist = "$bamlist.new";
        }
        
        my $newdir = $parsing->temp_filename;
        mkdir "$key_dir/$newdir" unless(-d "$key_dir/$newdir");
        chdir "$key_dir/$newdir";
        `split -l $line_num2split $bamlist`;
        chdir "$key_dir";
        opendir (CURR, "$key_dir/$newdir");
        my @fs = grep{/^x[a-z]{2}$/} readdir(CURR);
        closedir (CURR);
        
        foreach my $f(sort @fs)
        {
            print "Download starting for bams in the file $f!\n";
            #dwnld_wgs_or_rna(bamlist,ker directory,bam output directory)
            dwnld_wgs_or_rna("$key_dir/$newdir/$f","$key_dir","$bam_output_dir",$dwld_cmd,"$samtools");
           if ($action eq "all")
            {
                launch_RNA_Seq_pileup("$key_dir/$newdir/$f","$cds_sorted_path","$mpileup_outdir","$bam_output_dir","$samtools");
            }
        }
    }
    elsif ($action eq "mpileups")
    {
        launch_RNA_Seq_pileup("$bamlist","$cds_sorted_path","$mpileup_outdir","$bam_output_dir","$samtools");
    }
}

sub launch_RNA_Seq_pileup
{
    #launch_RNA_Seq_pileup("bamlist","cds_sorted directory path","rna_mpileups directory","bam output directory");          
    my @cmds;
    my $bam_file = shift; #file containing bams
    my $self = $bam_file and $bam_file = shift if ref $bam_file;
    my ($cds_dir,$mpileup_outdir,$bam_output_dir,$samtools) = @_; #path to sorted bed files, path to directory where mpileups will be stored, path to directory where bams are stored, samtools command
    
    $cds_dir =~ s/\/$//;
    $cds_dir = realpath($cds_dir);
    
    open (LPI,"$bam_file");

    while(my $r = <LPI>)
    {
        chomp($r);
        my @a = split("\t",$r);      
        my $bam = $bam_output_dir."/".$a[0]."/".$a[0].".bam";
       if (!(-e "$bam.bai"))
        {
            next;
        }
        else
        {
            opendir (DIR_FH,$cds_dir);
            my @cds = readdir DIR_FH;
            @cds = grep{/$a[1]/}@cds;
            closedir(DIR_FH);
            #Runs a foreach loop because there may be TCGA IDs that are the same but with different sample types.
            foreach my $cd(@cds)
            {
                chomp($cd);
                $cd =~ s/\/$//;
                my $out_id = "$a[0].$cd";
                $out_id =~ s/\.bed//;
                $cd = $cds_dir.'/'.$cd;#It will be used by samtool mpileup in combine with comment -l
             
                unless(-e $bam)
                {
                    print "No bam for $bam!\n";
                    next;
                }
                unless(-e $cd)
                {
                    print "No cd for $cd!\n";
                    next;
                }
               if (-e "$mpileup_outdir/$out_id")
                {
                    unless(-s "$mpileup_outdir/$out_id" > 0)
                    {
                        #do not repeat mpileup
                        push @cmds, "$samtools mpileup -d 8000 -l $cd $bam > $mpileup_outdir/$out_id";
                    }
                }
                else
                {
                    push @cmds, "$samtools mpileup -d 8000 -l $cd $bam > $mpileup_outdir/$out_id";
                }
            }
        }
    }
    close (LPI);
    #Do downloading and do mpileups for RNA-Seq Bams only, as WGS data will do MCE in sub;
    mce_map
    {
        chomp($_);
        print "Going to submit: $_\n";
        `$_`;          
    }@cmds;
  
    #Comment out this line if you want to keep the downloaded bams.
    Delete_Files_Recursively("$bam_output_dir","\.bam");
}

#Downloads rna or wgs bams and indexes them
sub dwnld_wgs_or_rna
{
    my $infile = shift; #file containing bams
    my $self = $infile and $infile = shift if ref $infile;
    my ($key_dir,$bam_output_dir,$dwld_cmd,$samtools) = @_; #path to gdc key directory, path to directory where bams will be downloaded, download command, samtools command
    
    open (DW,$infile);
    
    my $current_dir = getcwd();
    $bam_output_dir =~ s/\/$//;
    
    my $token = ImportAlldataIntoVar("$key_dir/gdc.key");
    my @cmds;
    my $cmd;
    while(my $r = <DW>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print "Working on $a[0] for TCGA sample $a[1]!\n";
        mkdir "$bam_output_dir/$a[0]" unless(-d "$bam_output_dir/$a[0]");
        #if a bam exists and there is no index file, curl or aria2c will continue downloading the bam otherwise it will download the bam as new
        if (-f "$bam_output_dir/$a[0]/$a[0].bam" and (!(-f "$bam_output_dir/$a[0]/$a[0].bam.bai")))
        {
            if ($dwld_cmd eq "curl")
            {
                $cmd = "$dwld_cmd -C - --header \'X-Auth-Token: $token\'  \'https://gdc-api.nci.nih.gov/legacy/data/$a[0]\'  -s --output \'$bam_output_dir/$a[0]/$a[0].bam\' ";
            }
            else
            {
                #aria2c is the download command if curl was not specified as the command in the command line when running script
                $cmd = "$dwld_cmd -s 16 -x 16 -c --dir $bam_output_dir/$a[0]/ -o $a[0].bam --header \'X-Auth-Token: $token\' https://gdc-api.nci.nih.gov/legacy/data/$a[0]";
            }
        }
        else
        {
            if ($dwld_cmd eq "curl")
            {
                $cmd = "curl --header \'X-Auth-Token: $token\'  \'https://gdc-api.nci.nih.gov/legacy/data/$a[0]\'  -s --output \'$bam_output_dir/$a[0]/$a[0].bam\' ";
            }
            else
            {
                #aria2c is the download command if curl was not specified as the command in the command line when running script
                $cmd = "$dwld_cmd -s 16 -x 16 --dir $bam_output_dir/$a[0]/ -o $a[0].bam --header \'X-Auth-Token: $token\' https://gdc-api.nci.nih.gov/legacy/data/$a[0]";
            }
        }
        
        unless(-f "$bam_output_dir/$a[0]/$a[0].bam.bai")
        {
            push @cmds, "$cmd;$samtools index $bam_output_dir/$a[0]/$a[0].bam";
        }        
    }
    close (DW);
    mce_map
    {
        my $cmd = $_;
        my($c1,$c2) = split(";",$cmd);
        my $bam = [split(" ",$cmd)]->[-1];
        print STDERR "Working on $bam!\n";
        
    	unless(-f "$bam.bai")
        {
            print "Going to run:\n $c1\n$c2\n";
            `$c1`;`$c2`;
        }
        else
        {
            print "Going to run:\n $c2\n";
            `$c2`;
        }
       
        print "$_ is done.\n";
    }@cmds;
}

sub download_files_from_gdc
{
    my $infile = shift; #file containing a list of Genotype or Copy number estimate files to download
    my $self = $infile and $infile = shift if ref $infile;
    my ($key_dir,$geno_cnv_dir,$dwld_cmd) = @_; #path to gdc key, path to directory where either Genotypes or Copy number estimate will be downloaded,array type (Genotypes or Copy number estimate), download command
    
    open (DW,"$infile");
    
    my $current_dir = getcwd();
    $key_dir =~ s/\/$//;
    my $token = ImportAlldataIntoVar("$key_dir/gdc.key");
    my @cmds;
    
    while(my $r = <DW>)
    {
        chomp($r);
        my @a = split("\t",$r);
        my $cmd;
        print "Working on $a[0] for TCGA sample $a[1]!\n";
        if (-f "$geno_cnv_dir/$a[0].$a[1]")
        {
            if ($dwld_cmd eq "curl")
            {
                $cmd = "$dwld_cmd -C - --header \'X-Auth-Token: $token\' \'https://gdc-api.nci.nih.gov/legacy/data/$a[0]\' -s --output \'$geno_cnv_dir/$a[0].$a[1]\' ";
            }
            else
            {
                #aria2c is the download command if curl was not specified as the command in the command line when running script
                $cmd = "$dwld_cmd -s 16 -x 16 -c --dir \'$geno_cnv_dir\' -o $a[0].$a[1] --header \'X-Auth-Token: $token\' https://gdc-api.nci.nih.gov/legacy/data/$a[0]";
            }
        }
        else
        {
            if ($dwld_cmd eq "curl")
            {
                $cmd = "$dwld_cmd  --header \'X-Auth-Token: $token\' \'https://gdc-api.nci.nih.gov/legacy/data/$a[0]\' -s --output \'$geno_cnv_dir/$a[0].$a[1]\' ";
            }
            else
            {
                #aria2c is the download command if curl was not specified as the command in the command line when running script
                $cmd = "$dwld_cmd -s 16 -x 16 --dir \'$geno_cnv_dir\' -o $a[0].$a[1] --header \'X-Auth-Token: $token\' https://gdc-api.nci.nih.gov/legacy/data/$a[0]";
            }
        }
        
        print "The following code is submitted: \n$cmd\n\n";
        push @cmds,$cmd;
    }
    close (DW);
    
    mce_map
    {
        `$_`;
        print "$_ is done.\n";
    }@cmds;
}

sub Delete_Files_Recursively
{
    my $current_folder = shift; #directory where files will be removed
    my $self = $current_folder and $current_folder = shift if ref $current_folder; 
    my $rgx = shift; #regular expression
    
    $current_folder =~ s/\/$//;
    $current_folder = realpath($current_folder) if $current_folder =~ /^\./;
    my @all = dir($current_folder);
    
    my @files2delete;

    foreach my $item(@all)
    { 	
        push(@files2delete,$item) if $item =~ /$rgx/i;
    }
    
    @files2delete = map {$_ = "$current_folder/$_"}@files2delete;
 
    #delete all the files in the array.
    map
    {
        my $file = $_;
       if (unlink("$file") == 1)
        {
            print "$file deleted successfully.\n";
        }
        else
        {
            print "$file was not deleted.\n";
        }
    }@files2delete;
}

##################Delete_Files_Recursively sub##################
sub dir
{
    my $current_folder = shift;
    my @all;

    chdir $current_folder;

    #Get the all files and folders in the given directory.
    my @both = glob("*");

    my @folders;
    foreach my $item(@both)
    {
	if (-d $item)
        {
            #Get all folders into another array - so that first the files will appear and then the folders.
	    push(@folders,$item);
	}
        else
        {
            #If it is a file just put it into the final array.
	    push(@all,$item);
	}
    }

    foreach my $this_folder(@folders)
    {
	#Add the directory name to the return list - comment the next line if you don't want this feature.
	push(@all,"$this_folder/");
        
	#Continue calling this function for all the folders
	my $full_path = "$current_folder/$this_folder";
        
	my @deep_items = dir($full_path); # :RECURSION:
	foreach my $item(@deep_items)
        {
	    push(@all,"$this_folder/$item");
	}
    }
    return @all;
}
##################end of Delete_Files_Recursively sub##################

sub ImportAlldataIntoVar
{
    $/ = undef;
    my $file = shift;
    open (FH, "$file");
    my $all = <FH>;
    close FH;
    $/ = "\n";
    return $all;
}

1;
