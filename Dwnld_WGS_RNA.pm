#!/usr/bin/perl -w
package TCGA_Lib::Dwnld_WGS_RNA;

use strict;
use warnings;
use FileHandle;
use LWP::Simple;
use FindBin qw($Bin);
use File::Copy qw(copy);
use MCE::Map;
use Cwd;
use Cwd 'realpath';
my $TCGA_Analysis_Dir = realpath("../");
use lib "$Bin/";
use Parsing_Routines;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

my $parsing = TCGA_Lib::Parsing_Routines->new;

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
    my $disease_abbr = shift;
    my $self = $disease_abbr and $disease_abbr = shift if ref $disease_abbr;
    my $exp_strat = shift;#Experimental Strategy
    my $data_type;
    my @diseases = split(" ",$disease_abbr);
    my %dis_info_hash;
    
    foreach my $disease_abbr(@diseases)
    {
        #Go to this web to change hex into ASCII;
        #https://en.wikipedia.org/wiki/ASCII;
        #The gdc html link should be: https://gdc-api.nci.nih.gov/legacy/files?;
        my $url;
        if($exp_strat =~ /Genotyping array/)
        {
           $data_type = shift;
           $url = "https://gdc-api.nci.nih.gov/legacy/files?filters=%7B'op':'and','content':%5B%7B'op':'in','content':%7B'field':'cases.project.project_id','value':%5B'TCGA-$disease_abbr'%5D%7D%7D,%7B'op':'in','content':%7B'field':'files.data_type','value':%5B'$data_type'%5D%7D%7D%5D%7D";
        }
        elsif($exp_strat =~ /RNA-Seq/ || $exp_strat =~ /WGS/)
        {
            $url = "https://gdc-api.nci.nih.gov/legacy/files?filters=%7B'op':'and','content':%5B%7B'op':'in','content':%7B'field':'files.data_type','value':%5B'Aligned%20reads'%5D%7D%7D,%7B'op':'in','content':%7B'field':'files.experimental_strategy','value':%5B'$exp_strat'%5D%7D%7D,%7B'op':'in','content':%7B'field':'cases.project.project_id','value':%5B'TCGA-$disease_abbr'%5D%7D%7D%5D%7D";
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
           print "Now we are querying all $disease_abbr cancer samples from TCGA project:\n\n",$url,"\n\n";
           if(defined $data_type)
           {
             `curl \'$url\' > \'$disease_abbr.$data_type.result.txt\'`;
           }
            else
            {
           `curl \'$url\' > \'$disease_abbr.result.txt\'`;
           }
    }
}

#Gets the metadata data of the file and places the UUID in a Payload.txt file.
sub metadata_collect
{
    my $input_file = shift;
    my $self = $input_file and $input_file = shift if ref $input_file;
    my $output = shift;
    open(my $in,$input_file) or die "failed to open $input_file: $!\n";
    open(MQ,">$output") or die "Can't open $output for output: $!\n";
    chomp(my @uuids = <$in>);
    close $in;
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
    close(MQ);
}

sub metadata_ids
{
   my $input_file = shift;
   my $self = $input_file and $input_file = shift if ref $input_file;
    my $output = shift;
    open(my $in,$input_file) or die "failed to open $input_file: $!\n";
    open(PAY,">$output") or die "Can't open $output for output: $!";
    chomp(my @uuids = <$in>);
    close $in;
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
    close($in);
    close(PAY);
}

#Uses the metadata file to get data(i.e. UUID, TCGA ID) and prints the data for the UUIDs to a datatable file.
sub parse_patient_id
{
    ################################################################
    # Parses the tumour-normal tissue number out of the Patient ID # 
    ################################################################
    
    my $input_file = shift;
    my $self = $input_file and $input_file = shift if ref $input_file;
    my $output_file = shift;
    open(my $in, $input_file) or die "failed to open $input_file: $!\n"; #open input_file
    chomp(my @data = <$in>);
    close $in;
    open(PID,">$output_file") or die "failed to open $output_file: $!\n";
    shift @data;  #remove headers
    
    foreach my $lines(@data)
    {
        my @split_line = split("\t",$lines);
        my $Patient_ID_column = $split_line[5];
        $Patient_ID_column = substr $Patient_ID_column, -3, 2;
        if($Patient_ID_column == 01)
        {
            $Patient_ID_column = 1;
        }
        if($Patient_ID_column == 02)
        {
            $Patient_ID_column = 2;
        }
        #.dtable(file.id | Patient ID | Sample Type | Tumour/Normal(num) | Platform)
        print PID $split_line[3],"\t", $split_line[1],"\t",$Patient_ID_column,"\t",$split_line[4],"\t",$split_line[0],"\n";
    }
    close(PID);
}

sub parse_meta_id
{
    my $meta_table=shift;
    my $self = $meta_table and $meta_table = shift if ref $meta_table;
    my $UUID_list = shift;
    my $meta_ids = shift; 
    my @ID_column;
    
    open(UUIDS, "$UUID_list") or die;
    chomp(my @all_bam_uuids = <UUIDS>);
    close UUIDS;
    
    open(my $fh,"$meta_table") or die "couldn't find file";
    chomp(my @table = <$fh>);
    close $fh;
    open(PMI,">$meta_ids") or die "Can't open file for output: $!";
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
                $bam_uuid = $id if($id eq $i);
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
    close(PMI);
}

sub ref_parse
{
    my %refs;
    my $meta_ID = shift;
    my $self = $meta_ID and $meta_ID = shift if ref $meta_ID;
    my $gdc_key_path = shift;
    my $reference = shift;
    my $token = ImportAlldataIntoVar("$gdc_key_path/gdc.key");
    
    open(my $fh,$meta_ID) or die "couldn't find file $meta_ID: $!\n";
    my @table = <$fh>;
    close($fh);
    open(RP,">$reference")or die "Can't open file $reference for output: $!\n";
    mkdir "ref_tmp";
   mce_map
   {
        chomp($_);
        my($refID,$UUID) = ($_ =~ /(\S*)\t(\S*)/);
        my $curl_cmd = "curl -# --header \'X-Auth-Token: $token\' \'https://gdc-api.nci.nih.gov/legacy/data/$refID\' --output \'ref_tmp/$refID.$UUID\' ";
        `$curl_cmd`;
    }@table;
   opendir REF, "ref_tmp" or die;
    my @ref_txts = grep{!/^\./ && -f "ref_tmp/$_"} readdir(REF);
    #@ref_txts = grep{!/^\./}@ref_txts;
    closedir(REF);
    
    foreach my $txt(@ref_txts)
    {
        my($RID,$UID) = split("\\.",$txt);
        $/ = undef;
        open(JSON,"ref_tmp/$txt") or die "$!";
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
    `rm -r ref_tmp`;
    foreach my $key(keys %refs)
    {
        print RP "$key\t$refs{$key}\n";
    }
    close(RP);
}

sub index_ids
{
    my $input_file = shift;
    my $self = $input_file and $input_file = shift if ref $input_file;
    my $output = shift;
    open(my $in,$input_file) or die "failed to open file $input_file: $!\n";
    chomp(my @uuids = <$in>);
    close $in;
    open(INDEX,">$output") or die "Can't open file $output for output: $!\n";
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
    close(INDEX);
}

sub pull_matched_tn_GDC
{
    my $input = shift;
    my $self = $input and $input = shift if ref $input;
    my $output = shift;
    
    my $total_tn = 0;
    open(MI,"$input") or die "Can't open file for input: $!\n";
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
    close(MI);    
    open(MO,">$output") or die "Can't open file $output for output: $!\n";
    my @all_keys = keys %ID_ref_platfm;
    foreach my $k(keys %uniq_ID_ref_platfm)
    {
        my @targets = grep{/$k/i;}@all_keys;
        if(@targets > 1)
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
                if($vv->[1] < 10)
                {
                    foreach my $xx(@tmp)
                    {
                        if($xx->[1] > 9)
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
    close(MO);
    if($total_tn == 0)
    {
        print STDERR "No tumor-normal pairs in your file: $input! \n";
    } 
}

sub Dwld_WGSBam_and_do_mpileup
{
    my $bamlist = shift;
    my $self = $bamlist and $bamlist = shift if ref $bamlist;
    my $key_dir = shift;
    my $wgs_output_dir = shift;
    my $ref_fullpath = shift;#./GRCh37-lite.fa
    my $mpileup_outdir = shift;#./mpileups
    my $action = shift;
    my $pairs = shift;
    my $VarScan_Path = shift;
    my $disease_abbr = shift;
    my $alt_ref = shift;
    my $tables = shift;
    
    $mpileup_outdir =~ s/\/$//;
    
    print STDERR "the bamlist dose not exist: $bamlist" and exit unless (-e "$bamlist");
    
    if($action eq "all" or $action eq "download")
    {
        if (-f ("$tables/already_done_WGS.txt"))
        {
            #Checks if there a NaN in the bamlist file and removes them as it will interfere with the code below.
            open(BAMI,"$bamlist") or die "Can't open $bamlist: $!\n";
            open(BAMO,">$bamlist.parse") or die "Can't open $bamlist: $!\n";
            
            while (my $r = <BAMI>)
            {
                chomp($r);
                
                 if ($r =~ /NaN/g)
                {
                    $r =~ s/NaN//;
                }
                print BAMO $r,"\n";   
            }
            close(BAMI);
            close(BAMO);
            
            #Matches the tumor bams in the bam list
            $parsing->vlookup("$bamlist.parse",1,"$tables/already_done_WGS.txt",1,1,"y","wgs_tum.txt");
            #greps all that didn't match
            `cat wgs_tum.txt|grep NaN > wgs_NaN.txt`;
            #Gets all columns except for the NaN column
            $parsing->pull_column("wgs_NaN.txt","1,2,3,4,5,6,7","wgs_get_norm.txt");
            #Matches the normal bams
            $parsing->vlookup("wgs_get_norm.txt",1,"$tables/already_done_WGS.txt",2,2,"y","wgs_norm.txt");
            #Gets the bams that did not meet any of the above criterea
            `cat wgs_norm.txt|grep NaN > $bamlist.new`;
            `rm wgs_*.txt`;
        }
        elsif(-f("$key_dir/already_done_WGS.txt"))
        {
            #Checks if there a NaN in the bamlist file and removes them as it will interfere with the code below.
            open(BAMI,"$bamlist") or die "Can't open $bamlist: $!\n";
            open(BAMO,">$bamlist.parse") or die "Can't open $bamlist: $!\n";
            
            while (my $r = <BAMI>)
            {
                chomp($r);
                
                 if ($r =~ /NaN/g)
                {
                    $r =~ s/NaN//;
                }
                print BAMO $r,"\n";   
            }
            close(BAMI);
            close(BAMO);
            
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
        mkdir "$newdir";
        opendir(CURR, "$newdir") or die "Can not open directory $newdir: $!";
        my @fs = readdir(CURR);
           @fs = grep{/^x[a-z]{2}$/}@fs;
        unlink @fs;
        closedir(CURR);
        chdir "$newdir" or die "Can not change into the dir $newdir: $!\n";
        
        #only split 1 pair into each file;
        `split -l 2 $bamlist.new`;
        chdir $key_dir;
        #Get all split files in the $newdir;
        opendir(CURR, "$newdir") or die "Can not open current dir: $!";
           @fs = readdir(CURR);
           @fs = grep{/^x[a-z]{2}$/}@fs;
        closedir(CURR);
        
        #Going to do downloading and do mpileups for these bams;
        my %PIDs;
        my $tag = 0;#record how many pairs submitted!
        
        foreach my $f(sort @fs)
        {
            #if there are more pairs than that of specified,
            #just wait it to finish;
            if ($tag > $pairs)
            {
                my $running_pids = $parsing->Get_PIDs();
                chomp($running_pids);
                my @rpids = split(" ",$running_pids);
                my @wgs_pids = keys %PIDs;
                my @wgs_left = Array_Intersect(\@rpids,\@wgs_pids);
                #wait here;
                while (@wgs_left >= $pairs)
                {
                    my $t = 3000;
                    print STDERR "Going to sleep for $t seconds\n";
                    sleep($t);
                    $running_pids = $parsing->Get_PIDs();
                    chomp($running_pids);
                    @rpids = split(" ",$running_pids);
                    @wgs_left = Array_Intersect(\@rpids,\@wgs_pids);
                }
            }
            my $PID_tmp;
            if ($action eq "all")
            {
                #Remove bams if their corresponding pid has gone;
                #Get running PIDs again;
                $PID_tmp = $parsing->Get_PIDs();
                chomp($PID_tmp);
                foreach my $p (keys %PIDs)
                {
                    unless ($PID_tmp=~/\b$p\b/)
                    {
                        #delete these bams
                        open my $L,"$PIDs{$p}" or die "Can not open the bamlist: $!";
                        #delete hash key $p and its value;
                        delete $PIDs{$p};
                        chomp(my @fs = <$L>);
                        map
                        {
                            my $f = $_;
                            $f =~ s/^([^\t]+).*/$1/;
                            print STDERR "PID $p has gone!\n",
                            "Going to delete bam file: $wgs_output_dir/$f/$f.* \n";
                            `rm -f $wgs_output_dir/$f/$f.*`;
                        }@fs;
                        close $L;   
                    }
                }
            }
            
            print "Now start to download bams in the file $f!\n";
            #dwnld_wgs_or_rna(bamlist,key directory,bam output directory)
            dwnld_wgs_or_rna("$key_dir/$newdir/$f","$key_dir","$wgs_output_dir");
            
            my $excluded_pids = $parsing->Get_PIDs();
            chomp($excluded_pids);
           
           if($action eq "all")
            {
                my @bam_pairs = mk_files_for_wgs("$key_dir/$newdir/$f","$wgs_output_dir");
                
                launch_wgs_mpileup_and_VarScan("$ref_fullpath","$mpileup_outdir",\@bam_pairs,$wgs_output_dir,$VarScan_Path,$disease_abbr,$alt_ref,"$tables");
            }
            
            my $npid = $parsing->Get_PIDs($excluded_pids);
            chomp($npid);#Need to remove newline!;
            if ($npid eq "")
            {
                print STDERR "No PID exists for the submitted job, maybe this computer can not submit the job!\n";
                exit;
            }
            
            $PIDs{$npid} = "$key_dir/$newdir/$f" if $npid =~ /^\d+$/;
            
            #Check the above submitted PID again;
            #In case of gdc.key expired
            #the above PID will be gone very quickly;
            #First get all PIDs again;
            print STDERR "Going to sleep 60s and check the submitted job $npid again\n" and sleep(60);
            my $TMP_PIDs = $parsing->Get_PIDs();
            chomp($TMP_PIDs);
            unless ($TMP_PIDs =~ /\b$npid\b/)
            {
                print STDERR "The gdc.key may be expired, as the PID is gone\n","Please check these bams included in the file: $newdir/$f\n";
                exit;                   
            }
            #record how many bam pairs have been submitted; 
            $tag++;
        }
    }
    elsif($action eq "mpileups")
    {
        my @bam_pairs = mk_files_for_wgs("$bamlist","$wgs_output_dir");
                
        launch_wgs_mpileup_and_VarScan("$ref_fullpath","$mpileup_outdir",\@bam_pairs,$wgs_output_dir,$VarScan_Path,$disease_abbr,$alt_ref,"$tables");
    }
}

sub Dwld_RNASeq_Bam_and_do_mpileup
{
    my $bamlist = shift;
    my $self = $bamlist and $bamlist = shift if ref $bamlist;
    my $key_dir = shift;#directory for gdc key and the temp directory for the file for downloading bams
    my $bam_output_dir = shift;
    my $ref_fullpath = shift;#./GRCh37-lite.fa
    my $mpileup_outdir = shift;#./rna_mpileups
    my $cds_sorted_path = shift;
    my $action = shift;
    my $line_num2split = shift;
    
    $bamlist =~ s/\/$//;
    print STDERR "bamlist $bamlist does not exist!\n" and exit unless(-e "$bamlist"); 
    $key_dir =~ s/\/$//;      
    print STDERR "Dir $key_dir does not exist!\n" and exit unless(-d "$key_dir");  
    $mpileup_outdir =~ s/\/$//;
    `mkdir -p $mpileup_outdir` unless(-d "$mpileup_outdir");
    
    if($action eq "all" or $action eq "download")
    {
        my $BAMs;
        if ($action eq "all")
        {
            #Filters through the table to find the matching bed files in cds_sorted and only includes those entries that have a matching bed in the cds_sorted directory.
            `ls $cds_sorted_path > cds_beds.txt`;
            open(CDS,"cds_beds.txt") or die "Can't open cds_beds.txt: $!\n";
            open(CDO,">cds_beds_out.txt") or die "Can't open cds_beds_out.txt: $!\n";
            while(my $r = <CDS>)
            {
                chomp($r);
                $r = substr($r,0,12);
                print CDO $r, "\n";
            }
            close(CDO);
            close(CDS);
            
            $parsing->vlookup("$bamlist",2,"cds_beds_out.txt",1,1,"y","cds_beds_lookup.txt");
            `grep -v NaN cds_beds_lookup.txt > cds_beds_new.txt`;
            my $bamlist_cds = "cds_beds_new.txt";
            open($BAMs,"$bamlist_cds") or die "can not open the bamlist file $bamlist: $!\n";
        }
        else
        {
            open($BAMs,"$bamlist") or die "can not open the bamlist file $bamlist: $!\n";
        }
        
        #Get already done pileup files and remove it from the $bamlist;
        print STDERR "Going to check already done pileups in the $mpileup_outdir\n";
        my @already_done = `ls $mpileup_outdir|grep '.'`;
           @already_done = grep{chomp;-s "$mpileup_outdir/$_";}@already_done;
        my %already_done = map{my($k,$v) = $_ =~ /([^\.]+)\.([^\.]+)/;$k=>$v}@already_done;
        open(my $NBAMs,">$bamlist.new") or die "can not write data into the file $bamlist.new: $!\n";
        while(my $b = <$BAMs>)
        {
            chomp($b);
            my @es = split("\t",$b);
            if(defined $already_done{$es[0]})
            {
                print STDERR "The TCGA ID $es[0].$es[1] has pileup result in the $mpileup_outdir and will not be submitted for pileup again\n";
            }
            else
            {
                print $NBAMs $b,"\n" ;          
            }
        }
        close $BAMs;
        close $NBAMs;
        `cp $bamlist $bamlist.backup -f`;
        $bamlist = "$bamlist.new";
        
        my $newdir = $parsing->temp_filename;
        mkdir "$key_dir/$newdir";
        opendir(CURR, "$key_dir/$newdir") or die "Can not open directory $key_dir/$newdir: $!";
        my @fs=readdir(CURR);
           @fs=grep{/^x[a-z]{2}$/}@fs;
        unlink @fs;
        closedir(CURR);
        chdir "$key_dir/$newdir";
        `split -l $line_num2split $bamlist`;
        chdir "../";
        opendir(CURR, "$key_dir/$newdir") or die "Can not open directory $key_dir/$newdir: $!";
           @fs = readdir(CURR);
           @fs = grep{/^x[a-z]{2}$/}@fs;
        closedir(CURR);
        
        foreach my $f(sort @fs)
        {
            print "Download starting for bams in the file $f!\n";
            #dwnld_wgs_or_rna(bamlist,ker directory,bam output directory)
            dwnld_wgs_or_rna("$key_dir/$newdir/$f","$key_dir","$bam_output_dir");
            if($action eq "all")
            {
                launch_RNA_Seq_pileup("$key_dir/$newdir/$f","$cds_sorted_path","$mpileup_outdir","$bam_output_dir");
            }
        }
    }
    elsif($action eq "mpileups")
    {
        launch_RNA_Seq_pileup("$bamlist","$cds_sorted_path","$mpileup_outdir","$bam_output_dir");
    }
}

sub launch_RNA_Seq_pileup
{
    #launch_RNA_Seq_pileup("bamlist","cds_sorted directory path","rna_mpileups directory","bam output directory");          
    my @cmds;
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $cds_dir = shift;
    my $mpileup_outdir = shift;
    my $bam_output_dir = shift;
    
    $cds_dir =~ s/\/$//;
    $cds_dir = realpath($cds_dir);
    
    open(LPI,"$infile") or die "Can't open $infile: $!\n";

    while(my $r = <LPI>)
    {
        chomp($r);
        my @a = split("\t",$r);      
        my $bam = $bam_output_dir."/".$a[0]."/".$a[0].".bam";
        if(!(-e "$bam.bai"))
        {
            next;
        }
        else
        {
            opendir(DIR_FH,$cds_dir) or die "Can't open directory $cds_dir: $!\n";
            my @cds = readdir DIR_FH;
            @cds = grep{/$a[1]/}@cds;
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
                    print STDERR "no bam for $bam: $!\n";
                    next;
                }
                unless(-e $cd)
                {
                    print STDERR "no cd for $cd\n";
                    next;
                }
                if(-e "$mpileup_outdir/$out_id")
                {
                    unless(-s "$mpileup_outdir/$out_id" > 0)
                    {
                        #do not repeat mpileup
                        push @cmds, "samtools mpileup -d 8000 -l $cd $bam > $mpileup_outdir/$out_id";
                    }
                }
                else
                {
                    push @cmds, "samtools mpileup -d 8000 -l $cd $bam > $mpileup_outdir/$out_id";
                }
            }
        }
    }
    close(LPI);
    #Do downloading and do mpileups for RNA-Seq Bams only, as WGS data will do MCE in sub;
    mce_map
    {
        chomp($_);
        print STDERR "Going to submit: $_\n";
        `$_`;          
    }@cmds;
  
    #Delete_Files_Recursively(searchingDir,file_rgx);
    #Comment out this line if you want to keep the downloaded bams.
    Delete_Files_Recursively("$bam_output_dir","\.bam");
}

sub Delete_Files_Recursively
{
    my $current_folder = shift;
    my $self = $current_folder and $current_folder = shift if ref $current_folder; 
    my $rgx = shift;
    
    $current_folder =~ s/\/$//;
    $current_folder = realpath($current_folder) if $current_folder =~ /^\./;
    my @all = dir($current_folder);
    
    my @files2delete;

    foreach my $item(@all)
    { 	
        push(@files2delete,$item) if $item =~ /$rgx/i;
    }
    
    @files2delete=map {$_ = "$current_folder/$_"}@files2delete;
 
    #delete all the files in the array.
    map
    {
        my $file = $_;
        if(unlink("$file") == 1)
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

    chdir($current_folder) or die("Cannot access folder $current_folder");

    #Get the all files and folders in the given directory.
    my @both = glob("*");

    my @folders;
    foreach my $item(@both)
    {
	if(-d $item)
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

#Downloads rna or wgs bams and indexes them
sub dwnld_wgs_or_rna
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $key_dir=shift;
    my $output_dir = shift;
    
    open(DW,$infile) or die "Can't open file for input: $!\n";
    
    my $current_dir = getcwd();
    $output_dir =~ s/\/$//;
    
    my $token = ImportAlldataIntoVar("$key_dir/gdc.key");
    my @cmds;
    my $cmd;
    while(my $r = <DW>)
    {
        chomp($r);
        my @a = split("\t",$r);
        print STDERR "working on $a[0] for TCGA sample $a[1]!\n";
        mkdir "$output_dir/$a[0]" unless(-d "$output_dir/$a[0]");
        #if a bam exists and there is no index file, curl will continue downloading the bam otherwise it will download the bam as new
        if (-f "$output_dir/$a[0]/$a[0].bam" and (!(-f "$output_dir/$a[0]/$a[0].bam.bai")))
        {
            $cmd = "curl -C - --header \'X-Auth-Token: $token\'  \'https://gdc-api.nci.nih.gov/legacy/data/$a[0]\'  -s --output \'$output_dir/$a[0]/$a[0].bam\' ";
        }
        else
        {
            $cmd = "curl --header \'X-Auth-Token: $token\'  \'https://gdc-api.nci.nih.gov/legacy/data/$a[0]\'  -s --output \'$output_dir/$a[0]/$a[0].bam\' ";
        }
        
        unless(-f "$output_dir/$a[0]/$a[0].bam.bai")
        {
            push @cmds, "$cmd;samtools index $output_dir/$a[0]/$a[0].bam";
        }        
    }
    close(DW);
    mce_map
    {
        my $cmd = $_;
        my($c1,$c2) = split(";",$cmd);
        my $bam = [split(" ",$cmd)]->[-1];
        print STDERR "working on $bam!\n";
        
    	if (-f "$bam")
        {
            print STDERR "Going to run:\n $c2\n";
            `$c2`;
    	}
        else
        {
            print STDERR "Going to run:\n $c1\n";
            `$c1`;
            print STDERR "Going to run:\n $c2\n";
            `$c2`;
    	}
        
        unless(-f "$bam.bai")
        {
            print STDERR "There is no index file for this bam. Going to run:\n $c1\n$c2\n";
            `$c1`;`$c2`;
        }
        print STDERR "$_ is done.\n";
    }@cmds;
   
   #Checks for index files of bams and if there are none for a scpecific bam it will go and redownload them as they are incomplete.
    my @bams = `ls $output_dir`;
    my @dwnld_bam;
    my $ctr = 0;
    for(my $i = 0;$i < scalar(@bams);$i++)
    {
        chomp($bams[$i]);
        opendir my $dir,"$output_dir/$bams[$i]";
        my @get_bai = readdir $dir;
        my @get_bam = grep{/\.bam$/}@get_bai;
        if(@get_bam)
        {
            @get_bai = grep{/\.bam.bai/}@get_bai;
            if(@get_bai)
            {
                next;
            }
            else
            {
                for(my $c = 0;$c < scalar(@cmds);$c++)
                {
                    my @a = split(" ",$cmds[$c]);
                    $a[9] = [split("/",$a[9])]->[-1];
                    if($get_bam[0] =~ /$a[9]/)
                    {
                        $dwnld_bam[$ctr] = $cmds[$c];
                        $ctr++;
                    }
                }
            }
        }
        else
        {
            next;
        }
    }
    #will only download bams if they were incomplete.
    if(@dwnld_bam)
    {
        mce_map
        {
            my $cmd = $_;
            my($c1,$c2) = split(";",$cmd);
            my $bam = [split(" ",$cmd)]->[-1];
            print STDERR "working on $bam!\n";
            if(-f "$bam")
            {
                print STDERR "Going to run:\n $c2\n";
                `$c2`;
            }
            else
            {
                print STDERR "Going to run:\n $c1\n";
                `$c1`;
                print STDERR "Going to run:\n $c2\n";
                `$c2`;
            }
            unless(-f "$bam.bai")
            {
                print STDERR "There is no index file for this bam. Going to run:\n $c1\n$c2\n";
                `$c1`;`$c2`;
            }
            print STDERR "$_ is done.\n";
            
        }@dwnld_bam;
    }
}

sub download_files_from_gdc
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $key_dir = shift;
    my $output_dir = shift;
    my $data_type = shift;
    
    open(DW,"$infile") or die "Can't open file for input: $!\n";
    
    my $current_dir = getcwd();
    $key_dir =~ s/\/$//;
    my $token = ImportAlldataIntoVar("$key_dir/gdc.key");
    my @cmds;
    
    while(my $r = <DW>)
    {
        chomp($r);
        my @a = split("\t",$r);
        my $cmd;
        print STDERR "working on $a[0] for TCGA sample $a[1]!\n";
        if (-f "$output_dir/$a[0].$a[1]")
        {
            $cmd = "curl -C - --header \'X-Auth-Token: $token\' \'https://gdc-api.nci.nih.gov/legacy/data/$a[0]\' -s --output \'$output_dir/$a[0].$a[1]\' ";
        }
        else
        {
            $cmd = "curl  --header \'X-Auth-Token: $token\' \'https://gdc-api.nci.nih.gov/legacy/data/$a[0]\' -s --output \'$output_dir/$a[0].$a[1]\' ";
        }
        
        print STDERR "The following code is submitted: \n$cmd\n\n";
        push @cmds,$cmd;
    }
    close(DW);
    
    mce_map
    {
        `$_`;
        print STDERR "$_ is done.\n";
    }@cmds;
}

sub mk_files_for_wgs
{
    my $infile = shift;
    my $self = $infile and $infile = shift if ref $infile;
    my $wgs_bam_dir = shift;
    
    open(FFW,$infile) or die "Can't open file for input: $!\n";
    
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
            print STDERR "not norm and tumor for the following line: \n", $n,"\n",
            $r, "\n";
            next;
        }      
        
        if (-e("$wgs_bam_dir/$tum/$tum.bam.bai") and -e("$wgs_bam_dir/$norm/$norm.bam.bai"))
        {
            my $T_bam="$wgs_bam_dir/$tum/$tum.bam";
            my $N_bam="$wgs_bam_dir/$norm/$norm.bam";
        
            if(-f "$T_bam" and -f "$N_bam")
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
    my $fasta = shift;
    my $self = $fasta and $fasta = shift if ref $fasta;
    my $mpileups_out_dir = shift;
    my $bam_fullpaths_ref = shift;
    my $wgs_output_dir = shift;
    my $VarScan_Path = shift;
    my $disease_abbr = shift;
    my $alt_ref = shift;
    my $tables = shift;
    
    if (-e "$tables/already_done_WGS.txt")
    {
        open(DON,">>$tables/already_done_WGS.txt") or die "Can't open file already_done_WGS.txt: $!\n";
    }
    else
    {
        open(DON,">>already_done_WGS.txt") or die "Can't open file already_done_WGS.txt: $!\n";
    }
    
    my @done_bams;
    my $tumor;
    my $normal;
    foreach my $tn_bam(@$bam_fullpaths_ref)
    {
        chomp($tn_bam);
        next if $tn_bam =~ /^\s*$/;
        my @a = split("\t",$tn_bam);
        
        $a[0] =~ /([^\/]+)\/[^\/]+\.bam$/;#revised!
        my $norm = $1;
        $norm .= ".".$a[2];
        print $norm," has been submitted\n";
        
        my $reference=$fasta;
        my $normal_bam=$a[0];
        my $tumor_bam=$a[1];
        
        #Check the size of bam;
        my $Gb=1024 * 1024;
        if((-s "$tumor_bam") < $Gb and (-s "$normal_bam") < $Gb)
        {
            print STDERR "Your bams are less than 1Gb, which is abnormal\n";
            print STDERR "tumor_bam: $tumor_bam\n";
            print STDERR "normal_bam: $normal_bam\n";
            exit;#stop all mpileups and need to check GDC, which may be under maintanence;
        }
        
        #Important!
        #Only with the right fasta ref for WGS Bams,
        #Varscan can give correct results;
        
        #Check compatibility between fasta ref and ref of WGS Bams;
        my $Ref_chr_label=`head -n 1 $reference`;
        print STDERR "fasta header is: \n $Ref_chr_label";
        my ($chr_info) = $Ref_chr_label =~ />((chr)*\d+)\s/msi;
        my $chr_tag;
        print STDERR "This is the parsed chr info \"$chr_info\" ","for $reference\n";
        if($chr_info =~ /chr/i)
        {
           $chr_tag = 0;
        }
        else
        {
           $chr_tag = 1;
        }
        
        my $tum_ref_tag = Bam_Ref_Checker($tumor_bam);
        print STDERR "tumor ref tag is $tum_ref_tag\n";
        my $norm_ref_tag = Bam_Ref_Checker($normal_bam);
        print STDERR "normal ref tag is $norm_ref_tag\n";
        
        if ($chr_tag == $tum_ref_tag and $chr_tag == $norm_ref_tag)
        {
            print STDERR "Your fasta ref $reference is compatable with these two bams: $tumor_bam and $normal_bam\n";
        }
        else
        {
            #try to use alternative fasta;
            print STDERR "You fasta ref $reference is NOT compatable with these two bams: $tumor_bam and $normal_bam\n";
            $reference = $alt_ref;
            print STDERR "We will try to use this fasta ref $reference, in which the chroms are labeled with 'chr'!\n";
        } 
        my $NT_pileup = "samtools mpileup -q 1 -B -f $reference $normal_bam $tumor_bam";
        my $output="$mpileups_out_dir/$norm";
        
        #Checks for existing files and goes to the next pair if the file exists. Delete files that may be incomplete.
        if (-e "$mpileups_out_dir/$norm.*.*")
        {
            next;
        }
        
        my $cmd="bash -c \"java -jar $VarScan_Path somatic <($NT_pileup) $output.varscan --mpileup 1\"";
        print STDERR "Running: \n$cmd\n\n";
        `$cmd`;
        
        $tumor = [split("/",$tumor_bam)]->[-1];
        $normal = [split("/",$normal_bam)]->[-1];
        $tumor =~ s/.bam//;
        $normal =~ s/.bam//;
        
        push @done_bams, $tumor,"\t",$normal,"\n";
        
        my $prev = " ";
        foreach my $bams(@done_bams)
        {
            if ($bams ne $prev)
            {
                print DON $bams;
                $prev = $bams;
            }
            else
            {
                next;
            }
        }
    }
    close(DON);
}

######################launch_wgs_mpileup_and_VarScan sub#######################################
sub Bam_Ref_Checker
{
    my $Bam = shift;
    my $self = $Bam and $Bam = shift if ref $Bam;
    my $tag;
    unless (-e "$Bam")
    {
        print STDERR "Please make sure the fullpath of your Bam is correct: $Bam $!";
        exit;
    }
    else
    {
        my $ls=`samtools view -H $Bam`;
        #print $ls,"\n";
        my @as=split("\n",$ls);
        foreach my $chr_ref(@as)
        {
            my $a = $chr_ref;
            ($tag) = $a =~ /\bSN:(\S+)\s/;
            if (defined $tag)
            {
                if($tag=~/^chr/i)
                {
                   print STDERR "Your Bam is aligned to a ref with chrom labeled as 'Chr*'\n",
                                "Please use the corresponding ref!\n";
                   print "0";
                }
                else
                {
                   print STDERR "Your Bam is aligned to a ref with chrom labeled in number\n",
                                "Please use the corresponding ref!\n";
                   print "1";                        
                }
                last;    
            }
        }
    }
    return $tag;
}
######################end of launch_wgs_mpileup_and_VarScan sub#######################################

sub Array_Intersect
{
    my %seen = ();
    my $array_one = shift;
    my @ones = Remove_duplicated_elements(@$array_one);
    my $array_two = shift;
    my @twos = Remove_duplicated_elements(@$array_two);
    my @total;
    push @total,@ones;
    push @total,@twos;
    my @intersect = grep{$seen{$_}++}@total;
    return @intersect;  
}

sub ImportAlldataIntoVar
{
    $/ = undef;
    my $file = shift;
    print STDERR "The gdc key file doesn't exist!\n",
    "please check the file: \n$file\n" unless(-e $file);
    open(FH, "$file") or die "Can't open file: $!";
    my $all = <FH>;
    close FH;
    $/ = "\n";
    return $all;
}

sub Remove_duplicated_elements
{
   my %seen = ();
   my @array = @_;
   #Remove duplicated elements in the array.
   @array = grep{!$seen{$_}++}@array;
   return @array;
 }

1;