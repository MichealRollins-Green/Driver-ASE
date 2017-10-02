#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$Bin/..";
use Parsing_Routines;
use Dwnld_WGS_RNA;
use File::Copy;
use Getopt::Long;
use autodie;
use Cwd 'realpath';
use File::Basename;
no warnings 'once';

my $time = localtime;
print "Script started: $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;
#Creates objects of the modules that are going to be used in this script.
my $parsing = Driver_ASE_Lib::Parsing_Routines->new;
my $dwnld = Driver_ASE_Lib::Dwnld_WGS_RNA->new;

#the different oprions for this script. The out_dir and ase_dir are optional.
GetOptions(
    'cancer|c=s' => \my $cancer_type, #e.g. OV
    'Expstrategy|E=s' => \my $Exp_Strategy, #e.g. WGS RNA-Seq
    'option|o=s' => \my $option, #either all, download or mpileups
    'intersect|i=s' => \my $overlap, #choice to overlap IDs
    'number|n=i' => \my $number, #the number of bams to download
    'download|d=s' => \my $dwnld_cmd, #curl or aria2c (if aria or aria2 is entered, it changes them to aria2c as that is the command)
    'samtools|s=s' => \my $samtools, #path to samtools
    'VarScan|V=s' => \my $VarScan_Path, #path to the VarScan jar file
    'key|k=s'=> \my $key, #path to the gdc key file
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("3.0");

#If -h was specified on the command line then the usage of the program will be printed.
if ($help)
{
    $parsing->usage("3.0");
}

#If the disease abbr or file type was not specified then an error will be printed and the usage of the program will be shown.
if (!defined $cancer_type || !defined $Exp_Strategy)
{
    print STDERR "Cancer type and/or experimental strategy was not entered!\n";
    $parsing->usage("3.0");
}

if ($Exp_Strategy ne "RNA-Seq" and $Exp_Strategy ne "WGS")
{
    print STDERR "The experimental strategy entered is incorrect for this pipeline!\n";
    $parsing->usage("3.0");
}

if (!defined $samtools)
{
    $samtools = "samtools";
}

my $Driver_ASE_Dir = realpath("../../");
my $rna_dwnlds = "rna_dwnlds";
my $wgs_dwnlds = "wgs_dwnlds";
my $tables = "$cancer_type\_tables";
my $cds_sorted = "cds_sorted";
my $rna_mpileups = "rna_mpileups";
my $wgs_mpileups = "wgs_mpileups";
my $database_path = "$Driver_ASE_Dir/Database";
#Directory where all analysis data will be going in.
my $Analysispath = realpath("../../Analysis");
my $RNA_Path = "$Analysispath/$cancer_type/RNA_Seq_Analysis";
my $ase = "$RNA_Path/ase";
my $Intersect = "$cancer_type\_RNA_WGS_Geno.txt";
my $downloadtable = "final_downloadtable_$cancer_type";
my $RNA_table = "$downloadtable\_RNA-Seq.txt";
my $WGS_table = "$downloadtable\_WGS.txt";
my $Genotypes_table = "$cancer_type.Genotypes.id2uuid.txt";
my $reference_file = "$cancer_type\_reference.txt";
my $GRCh37_file = "GRCh37-lite.fa";
my $hg19_file = "hg19.fa";

#If the path to the gdc key was not specified then an error will be printed and the usage of the program will be shown.
if (!defined $key)
{
    print STDERR "GDC key fullpath was not entered or the fullpath to it was not correct!\n";
    $parsing->usage("3.0");
}

$parsing->check_directory_existence("$database_path","$key"); #check if directories or files exist
$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

if ($Exp_Strategy eq "WGS")
{
    if (!defined $VarScan_Path)
    {
	print STDERR "The VarScan path was not entered!\n";
	$parsing->usage("3.0");
    }
    else
    {
	$parsing->CheckSoftware("$VarScan_Path");
    }
}

#The default choice will be to do download and mpileup if no option was specified
#mpileups will only work if there are bams that have been downloaded and indexed
$option = "all" unless defined $option;
if ($option ne "all" and $option ne "download" and $option ne "mpileups")
{
    print STDERR "The choice must be either all, download or mpileups!\n";
    $parsing->usage("3.0");
}

#Defaults to curl if no download command was specified
if (!defined $dwnld_cmd)
{
    $dwnld_cmd = "curl";
    print "No download command specified! defaulting to $dwnld_cmd.\n";
}
elsif ($dwnld_cmd eq "curl" or $dwnld_cmd eq "aria2c")
{
    print "Using $dwnld_cmd as the download command.\n";
}
elsif ($dwnld_cmd eq "aria" or $dwnld_cmd eq "aria2")
{
    print "$dwnld_cmd entered, converting to aria2c.\n";
    $dwnld_cmd = "aria2c";
}
else
{
    print STDERR "The download command must be either curl or aria2c!\n";
    $parsing->usage("3.0")
}

$parsing->CheckSoftware("$samtools","$dwnld_cmd");

mkdir "$Driver_ASE_Dir/Analysis" unless (-d "$Driver_ASE_Dir/Analysis");

if ($Exp_Strategy eq "RNA-Seq")
{
    $option = $parsing->determine_option("$RNA_Path","$cds_sorted","$option");
}

if ($Exp_Strategy eq "RNA-Seq")
{
    $number = 45 unless defined $number;
}
elsif ($Exp_Strategy eq "WGS")
{
    if (!defined $number)
    {
        $number = 2;
    }
    elsif (0 != $number % 2)
    {
        print STDERR "$Exp_Strategy are downloaded in pairs! Enter an even number for the -n option";
        exit;
    }
}

if (!defined $overlap)
{
    $overlap = "no";
}
elsif (lc $overlap ne "n" and lc $overlap ne "no" and lc $overlap ne "y" and lc $overlap ne "yes")
{
    print STDERR "Incorrect option entered! options are n|no|y|yes.\n";
    print $parsing->usage("3.0");
}

#Makes the directory where all of the results will be processed and stored.
`mkdir -p $Analysispath/$cancer_type` unless (-d "$Analysispath/$cancer_type");
chdir "$Analysispath/$cancer_type";

my $OUT_DIR;
#defaults to a directory if no output directory was specified in the command line.

if ($Exp_Strategy eq "RNA-Seq")
{
    $OUT_DIR = "$Analysispath/$cancer_type/$rna_dwnlds/bams";
}
elsif ($Exp_Strategy eq "WGS")
{
    $OUT_DIR = "$Analysispath/$cancer_type/$wgs_dwnlds/bams";
}
else
{
    print STDERR "Experimental strategy is incorrect for this pipeline!\n";
    $parsing->usage("3.0");
}

`mkdir -p "$OUT_DIR"`;

$OUT_DIR =~ s/\/$//;
`mkdir -p "$OUT_DIR"` unless (-d "$OUT_DIR");
    
#If the entered was mpileups then it will check if the directory where the bams were downloaded exists.
if ($option eq "mpileups")
{   
    if ($Exp_Strategy eq "WGS" or $Exp_Strategy eq "RNA-Seq")
    {
        if (!(-d "$OUT_DIR"))
        {
            print STDERR "$OUT_DIR does not exist! No downloads will be performed.\n";
            print STDERR "Can't run mpileups on files that do not exits!\n";
            exit;
        }
    }  
}

my $rna_wgs_dir = dirname($OUT_DIR);
chdir "$rna_wgs_dir";

#Gets the full path to the directory if the path was not set to the default directory.
unless ($OUT_DIR eq "$Analysispath/$cancer_type/$rna_dwnlds/bams" or $OUT_DIR eq "$Analysispath/$cancer_type/$wgs_dwnlds/bams")
{
    $rna_wgs_dir = realpath("../$rna_wgs_dir");
    $OUT_DIR = realpath("../$OUT_DIR");
}

#Check gdc key and mv it to db first!
#copyfile2newfullpath(path to gdc key file,path where gdc key file will be copied)
$parsing->copyfile2newfullpath("$key","$rna_wgs_dir/gdc.key");

#Checks if a table file does not exist in the tables directory of the cancer type.
if (($Exp_Strategy eq "RNA-Seq" and !(-f "$Analysispath/$cancer_type/$tables/$RNA_table")) or ($Exp_Strategy eq "WGS" and !(-f "$Analysispath/$cancer_type/$tables/$WGS_table")))
{   
    #parses the gdc website manifest of the specified cancer type and prints to a results file.
    #gdc_parser($cancer_type(e.g. OV),$Exp_Strategy (RNA-Seq or WGS))
    $dwnld->gdc_parser($cancer_type,$Exp_Strategy);
    
    #Gets the metadata data of the file and places the UUID in a Payload.txt file.
    #metadata_collect(.result.txt file from gdc_parser,output file)
    $dwnld->metadata_collect("$cancer_type.result.txt","Payload.txt");
    
    #Gets metadata files for each UUID and prints to a metadata file.
    `curl --request POST --header "Content-Type: application/json" --data \@Payload.txt 'https://gdc-api.nci.nih.gov/legacy/files' > $cancer_type\_metadata.txt`;
    
    #Uses the metadata file to get data(i.e. UUID, TCGA ID) an prints the data for the UUIDs to a datatable file.
    #parse_patient_id(_metadata.txt file created from curl,output file)
    $dwnld->parse_patient_id("$cancer_type\_metadata.txt","$cancer_type.datatable.txt");
    
    #metadata_ids(.result.txt created from gdc_parser,output file)
    $dwnld->metadata_ids("$cancer_type.result.txt","Payload.txt");
    
    `curl --request POST --header "Content-Type: application/json" --data \@Payload.txt 'https://gdc-api.nci.nih.gov/legacy/files' > $cancer_type\_metadata.txt`;
    
    open (my $meta,"$cancer_type\_metadata.txt");
    open (ME,">$cancer_type.edit.metadata.txt");
    chomp(my @metaedit = <$meta>);
    close ($meta);
    for(my $i = 0;$i < scalar(@metaedit);$i++)
    {
        $metaedit[$i] =~ s/\r//g;
    }
    @metaedit = grep{!/\t\s+\t/}@metaedit;
    foreach(@metaedit)
    {
        print ME "$_\n";
    }
    close (ME);
    
    #pulls the UUIDs from the edit metadata file
    #pull_column(.edit.metadata.txt file,column(s) to pull,output file)
    $parsing->pull_column("$cancer_type.edit.metadata.txt","8","temp.UUID");
    
    #Strips the headers in the file.
    #strip_head(file to strip headers,output file)
    $parsing->strip_head("temp.UUID","$cancer_type.UUID.txt");
    `rm temp.UUID`;
    
    #parse_meta_id(.edit.metadata.txt file,output file from strip_head,output file)
    $dwnld->parse_meta_id("$cancer_type.edit.metadata.txt","$cancer_type.UUID.txt","meta_ids.txt");
    
    #ref_parse(output file from parse_meta_id,$rna_wgs_dir (directory where gdc key is located),output file)
    $dwnld->ref_parse("meta_ids.txt","$rna_wgs_dir","$reference_file",$dwnld_cmd);
    
    #index_ids(.result.txt file from gdc_parser,output file)
    $dwnld->index_ids("$cancer_type.result.txt","Payload.txt");
    
    `curl --request POST --header "Content-Type: application/json" --data \@Payload.txt 'https://gdc-api.nci.nih.gov/legacy/files' > index_file_ids.txt`;
    
    #vlookup(lookupFile,queryCol,sourceFile,lookupCol,returnCol(s),append(y/n),outputFile)
    
    #e.g. vlookup(lookupfile,3,sourcefile,4,"1,2,4,6","y",outputFile)

    #Will search each column 3 entry of lookupfile within column 4 of sourceFile and will append columns 1,2,4,6 of sourceFile to the end of each row in lookupfile.

    #N.B. only works on tab-delimited files
    $parsing->vlookup("$cancer_type.datatable.txt","1","$reference_file","1","2","y","temp");
    
    if ($Exp_Strategy eq "RNA-Seq")
    {
        $parsing->vlookup("temp","1","index_file_ids.txt","3","2","y","$downloadtable.txt");
    }
    else
    {
        $parsing->vlookup("temp","1","index_file_ids.txt","3","1","y","$downloadtable.txt");
    }
    
    `mkdir $Analysispath/$cancer_type/$tables` unless (-d "$Analysispath/$cancer_type/$tables");
 
    if ($Exp_Strategy eq "WGS")
    {
        `sort -k2,2 -k3,3 -k6,6 $downloadtable.txt > $downloadtable\_sorted.txt`;
        #No header in the output, thus no need to strip head!
        #pull_matched_tn_GDC($downloadtable_sorted.txt file,output file)
        $dwnld->pull_matched_tn_GDC("$downloadtable\_sorted.txt","$downloadtable\_$Exp_Strategy.txt");
	$parsing->vlookup("$downloadtable\_$Exp_Strategy.txt",1,"$cancer_type.result.txt",1,4,"y","$downloadtable\_$Exp_Strategy\_size.txt");
        open (SIZE,"$downloadtable\_$Exp_Strategy\_size.txt");
        open (SO,">$downloadtable\_$Exp_Strategy\_convert.txt");
        #convert the size of the WGS bams to gigabytes.
        while (my $r = <SIZE>)
        {
            chomp($r);
            my @con = split("\t",$r);
            my $vert = pop @con;
            $vert = $vert/1000/1000/1000;
            $vert = eval sprintf('%.2f',$vert);
            push(@con,$vert);
            my $size_wgs = join("\t",@con);
	    $size_wgs =~ s/\r//;
            print SO $size_wgs,"\n";
        }
        close (SIZE);
        close (SO);
	
        #filter table and remove bams aligning to NCBI36 or HG18;
        `cat $downloadtable\_$Exp_Strategy\_convert.txt | grep NCBI36 -v | grep -v HG18 > WGS_tmp.txt;mv WGS_tmp.txt $downloadtable\_$Exp_Strategy.txt`;
	copy("$downloadtable\_$Exp_Strategy.txt","$Analysispath/$cancer_type/$tables");
    }
    elsif ($Exp_Strategy eq "RNA-Seq")
    {
        #filter table and remove bams that may be aligned to HG18, HG19_Broad_variant, GRChr37-lite, NCBI36 or non aligned bams;
        `cat $downloadtable.txt|grep NCBI36 -v|grep -v HG18 |grep -v HG19_Broad_variant |grep -v GRCh37-lite > $downloadtable\_$Exp_Strategy.txt`;
        copy("$rna_wgs_dir/$downloadtable\_$Exp_Strategy.txt","$Analysispath/$cancer_type/$tables");
    }
}
else
{
    if ($Exp_Strategy eq "WGS")
    {
        print "It seems that there is a table in the directory $Analysispath/$cancer_type/$tables. $WGS_table\n";
    }
    else
    {
        print "It seems that there is a table in the directory $Analysispath/$cancer_type/$tables. $RNA_table\n";
    }
}

if (lc $overlap eq "y" || lc $overlap eq "yes")
{
    $parsing->check_directory_existence("$Analysispath/$cancer_type/$tables/$RNA_table","$Analysispath/$cancer_type/$tables/$WGS_table","$Analysispath/$cancer_type/$tables/$Genotypes_table"); #checks if RNA-Seq, WGS and Genotypes tables exist before overlapping them
    if (!(-s "$Analysispath/$cancer_type/$tables/$RNA_table" == 0) and !(-s "$Analysispath/$cancer_type/$tables/$WGS_table" == 0) and !(-s "$Analysispath/$cancer_type/$tables/$Genotypes_table" == 0))
    {
        $parsing->Overlap_RNA_WGS_Geno("$Analysispath/$cancer_type/$tables","$RNA_table","$WGS_table","$Genotypes_table","$Intersect","$Exp_Strategy","$cancer_type");
    }
    else
    {
        print "Either $RNA_table, $WGS_table and/or $Genotypes_table is empty!\n";
        exit;
    }
}

 if ($Exp_Strategy eq "WGS")
    {
        #Dwld_WGSBam_and_do_mpileup($Analysispath/$cancer_type/$tables/$downloadtable\_$Exp_Strategy.txt,$rna_wgs_dir (path to gdc key),$OUT_DIR (path to output directory),$database_path/$GRCh37_file (path to GRCh37-lite fa file),$rna_wgs_dir/$wgs_mpileups (path to directory where wgs mpileups will be stored),$option (all,download or mpileups),$number (pairs to split),$VarScan_Path (path to VarScan jar file),$cancer_type (e.g. OV),$database_path/$hg19_file (path to hg19 fa file),$Analysispath/$cancer_type/$tables (path to the table directory of the cancer type),$dwnld_cmd (download command (curl aria2c)),$samtools (path/command for samtools))
        mkdir "$rna_wgs_dir/$wgs_mpileups" unless (-d "$rna_wgs_dir/$wgs_mpileups");
        $dwnld->Dwld_WGSBam_and_do_mpileup("$Analysispath/$cancer_type/$tables/$downloadtable\_$Exp_Strategy.txt","$rna_wgs_dir","$OUT_DIR","$database_path/$GRCh37_file","$rna_wgs_dir/$wgs_mpileups",$option,$number,$VarScan_Path,$cancer_type,"$database_path/$hg19_file","$Analysispath/$cancer_type/$tables",$dwnld_cmd,"$samtools");
        copy("already_done_WGS.txt","$Analysispath/$cancer_type/$tables");
    }
    elsif ($Exp_Strategy eq "RNA-Seq")
    {
	if ($option ne "download")
	{
	    `mkdir -p $ase/$rna_mpileups` unless (-d "$ase/$rna_mpileups");
	}
	#Dwld_RNASeq_Bam_and_do_mpileup($Analysispath/$cancer_type/$tables/$downloadtable\_$Exp_Strategy.txt,$rna_wgs_dir (path to gdc key),$OUT_DIR (path to output directory),$database_path/$GRCh37_file (path to GRCh37-lite fa file),$ase/$rna_mpileups (path to the directory where rna mpileups will be stored),$cds_sorted (path to directory where cds sorted files are stored),$option (all,download or mpileups),$option (all,download or mpileups),$number(number of base to split into files),$dwnld_cmd (curl or aria2c),$samtools (path/command for samtools))
	$dwnld->Dwld_RNASeq_Bam_and_do_mpileup("$Analysispath/$cancer_type/$tables/$downloadtable\_$Exp_Strategy.txt","$rna_wgs_dir","$OUT_DIR","$database_path/$GRCh37_file","$ase/$rna_mpileups","$RNA_Path/$cds_sorted",$option,$number,$dwnld_cmd,"$samtools"); 
    }

chdir "$rna_wgs_dir";

my @del_files = $parsing->get_only_files_in_dir("$rna_wgs_dir");
#Removes all files after the processes have finished.
for(my $i = 0;$i < scalar(@del_files);$i++)
{
    `rm "$del_files[$i]"`;
}

`rm -rf tmp*`;

print "All jobs have finished for $cancer_type.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
