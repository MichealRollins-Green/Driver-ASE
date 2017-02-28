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

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;
#Creates objects of the modules that are going to be used in this script.
my $parsing = TCGA_Lib::Parsing_Routines->new;
my $dwnld = TCGA_Lib::Dwnld_WGS_RNA->new;

#the different oprions for this script. The out_dir and ase_dir are optional.
GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'exp_strat|e=s' => \my $Exp_Strategy,#e.g. WGS RNA-Seq
    'choice|c=s' => \my $choice,#either all, download or mpileups
    'number|n=i' => \my $number,#the number of bams to download
    'key|k=s'=> \my $key,
    'var_path|v=s' => \my $VarScan_Path,
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("3.0");

#If -h was specified on the command line then the usage of the program will be printed.
if($help)
{
    $parsing->usage("3.0");
}

my $TCGA_Pipeline_Dir = realpath("../../");
my $database_path = "$TCGA_Pipeline_Dir/Database";
#Directory where all analysis data will be going in.
mkdir "$TCGA_Pipeline_Dir/Analysis" unless(-d "$TCGA_Pipeline_Dir/Analysis");
my $Analysispath = realpath("../../Analysis");
my $RNA_Path = "$Analysispath/$disease_abbr/RNA_Seq_Analysis";
my $rna_dwnlds = "rna_dwnlds";
my $wgs_dwnlds = "wgs_dwnlds";
my $ase = "$RNA_Path/ase";
my $tables = "$disease_abbr\_tables";
my $cds_sorted = "cds_sorted";
my $rna_mpileups = "rna_mpileups";

#If the disease abbr or file type was not specified then an error will be printed and the usage of the program will be shown.
if(!defined $disease_abbr || !defined $Exp_Strategy)
{
    print "disease type and/or experimental strategy was not entered!\n";
    $parsing->usage("3.0");
}

#If the path to the gdc key was not specified then an error will be printed and the usage of the program will be shown.
if(!defined $key or (!(-f $key)))
{
    print "gdc key fullpath was not entered or the fullpath to it was not correct!\n";
    print $key,"\n";
    $parsing->usage("3.0");
}

if ($Exp_Strategy ne "RNA-Seq" and $Exp_Strategy ne "WGS")
{
    print STDERR "The the experimental strategy entered must be RNA-Seq or WGS as these are what this pipeline deals with\n";
    exit;
}

if($Exp_Strategy eq "WGS")
{
	if(!defined $VarScan_Path)
	{
		print STDERR "The VarScan path was not entered in\n";
		$parsing->usage("3.0");
	}
	elsif(!(-e $VarScan_Path))
	{
		print STDERR "$VarScan_Path does not exist, enter in the correct path\n";
	}
}

#The default choice will be to do download and mpileup if no option was specified
#mpileups will only work if there are bams that have been downloaded and indexed
$choice = "all" unless defined $choice;
if($choice ne "all" and $choice ne "download" and $choice ne "mpileups")
{
    print "The choice must be either all, download or mpileups.\n";
    $parsing->usage("3.0");
}

#Check if there is no Database directory
if(!(-d "$database_path"))
{
    print STDERR "$database_path does not exist, it was either moved, renamed, deleted or has not been downloaded.\nPlease check the README.md file on the github page to find out where to get the Database directory.\n";
    exit;
}

if ($Exp_Strategy eq "RNA-Seq")
{   
    if(-d $RNA_Path)
    {
        if (!(-d "$RNA_Path/cds_sorted"))
        {
            print STDERR "$RNA_Path/cds_sorted does not exist, only downloads can be performed!\n";
            exit;
        } 
    }
    else
    {
        print STDERR "$RNA_Path does not exist, only downloads will be performed!\n";
        $choice = "download";
    }
}

#Makes the directory where all of the results will be processed and stored.
`mkdir -p $Analysispath/$disease_abbr` unless(-d "$Analysispath/$disease_abbr");
chdir "$Analysispath/$disease_abbr";

if($Exp_Strategy eq "RNA-Seq")
{
    $number = 45 unless defined $number;
}
elsif($Exp_Strategy eq "WGS")
{
            $number = 2 unless defined $number;
}

my $OUT_DIR;
#defaults to a directory if no output directory was specified in the command line.

if($Exp_Strategy eq "RNA-Seq")
{
    $OUT_DIR = "$Analysispath/$disease_abbr/$rna_dwnlds/bams";
}
elsif($Exp_Strategy eq "WGS")
{
    $OUT_DIR = "$Analysispath/$disease_abbr/$wgs_dwnlds/bams";
}
else
{
    print STDERR "File type must be either RNA-Seq or WGS.\n";
    exit;
}

`mkdir -p "$OUT_DIR"`;

$OUT_DIR =~ s/\/$//;
`mkdir -p "$OUT_DIR"` unless(-d "$OUT_DIR");
    
#If the entered was mpileups then it will check if the directory where the bams were downloaded exists.
if($choice eq "mpileups")
{   
    if($Exp_Strategy eq "WGS" or $Exp_Strategy eq "RNA-Seq")
    {
        if(!(-d "$OUT_DIR"))
        {
            print "$OUT_DIR does not exist meaning that there were no downloads performed.\n";
            print "Can't run mpileups on files that do not exits.\n";
            exit;
        }
    }  
}

my $rna_wgs_dir = dirname($OUT_DIR);
chdir "$rna_wgs_dir" or die "Can't change to directory $rna_wgs_dir: $!\n";


#Gets the full path to the directory if the path was not set to the default directory.
unless($OUT_DIR eq "$Analysispath/$disease_abbr/$rna_dwnlds/bams" or $OUT_DIR eq "$Analysispath/$disease_abbr/$wgs_dwnlds/bams")
{
    $rna_wgs_dir = realpath("../$rna_wgs_dir");
    $OUT_DIR = realpath("../$OUT_DIR"); 
}

#Checks if a table file does not exist in the tables directory of the cancer type.
if(!(-f "$Analysispath/$disease_abbr/$tables/final_downloadtable_$disease_abbr\_$Exp_Strategy.txt"))
{
    #Check gdc key and mv it to db first!
    #copyfile2newfullpath(path to gdc key file,path where gdc key file will be copied)
    $parsing->copyfile2newfullpath("$key","$rna_wgs_dir/gdc.key");
    
    #parses the gdc website manifest of the specified cancer type and prints to a results file.
    #gdc_parser(cancer type(e.g. OV),type of data(RNA-Seq or WGS))
    $dwnld->gdc_parser($disease_abbr,$Exp_Strategy);
    
    #Gets the metadata data of the file and places the UUID in a Payload.txt file.
    #metadata_collect(.result.txt file from gdc_parser,output file)
    $dwnld->metadata_collect("$disease_abbr.result.txt","Payload.txt");
    
    #Gets metadata files for each UUID and prints to a metadata file.
    `curl --request POST --header "Content-Type: application/json" --data \@Payload.txt 'https://gdc-api.nci.nih.gov/legacy/files' > $disease_abbr\_metadata.txt`;
    
    #Uses the metadata file to get data(i.e. UUID, TCGA ID) an prints the data for the UUIDs to a datatable file.
    #parse_patient_id(_metadata.txt file created from curl,output file)
    $dwnld->parse_patient_id("$disease_abbr\_metadata.txt","$disease_abbr.datatable.txt");
    
    #metadata_ids(.result.txt,output file)
    $dwnld->metadata_ids("$disease_abbr.result.txt","Payload.txt");
    
    `curl --request POST --header "Content-Type: application/json" --data \@Payload.txt 'https://gdc-api.nci.nih.gov/legacy/files' > $disease_abbr\_metadata.txt`;
    
    open(my $meta,"$disease_abbr\_metadata.txt") or die "Can't open file $disease_abbr\_metadata.txt for input: $!";
    open(ME,">$disease_abbr.edit.metadata.txt") or die "Can't open file $disease_abbr.edit.metadata.txt for output: $!";
    chomp(my @metaedit = <$meta>);
    close($meta);
    for(my $i = 0;$i < scalar(@metaedit);$i++)
    {
        $metaedit[$i] =~ s/\r//g;
    }
    @metaedit = grep{!/\t\s+\t/}@metaedit;
    foreach(@metaedit)
    {
        print ME "$_\n";
    }
    close(ME);
    
    #pulls the UUIDs from the edit metadata file
    #pull_column(.edit.metadata.txt,column(s) to pull,output file)
    $parsing->pull_column("$disease_abbr.edit.metadata.txt","2","temp.UUID");
    
    #Strips the headers in the file.
    #strip_head(file to strip headers,output file)
    $parsing->strip_head("temp.UUID","$disease_abbr.UUID.txt");
    `rm temp.UUID`;
    
    #parse_meta_id(.edit.metadata.txt,output file from strip_head,output file)
    $dwnld->parse_meta_id("$disease_abbr.edit.metadata.txt","$disease_abbr.UUID.txt","meta_ids.txt");
    
    #ref_parse(output file from parse_meta_id,key directory,output file)
    $dwnld->ref_parse("meta_ids.txt","$rna_wgs_dir","reference.txt");
    
    #index_ids(.result.txt,output file)
    $dwnld->index_ids("$disease_abbr.result.txt","Payload.txt");
    
    `curl --request POST --header "Content-Type: application/json" --data \@Payload.txt 'https://gdc-api.nci.nih.gov/legacy/files' > index_file_ids.txt`;
    
    #vlookup(lookupFile,queryCol,sourceFile,lookupCol,returnCol(s),append(y/n),outputFile)
    
    #e.g. vlookup(lookupfile,3,sourcefile,4,"1,2,4,6","y",outputFile)

    #Will search each column 3 entry of lookupfile within column 4 of sourceFile and will append columns 1,2,4,6 of sourceFile to the end of each row in lookupfile.

    #N.B. only works on tab-delimited files
    $parsing->vlookup("$disease_abbr.datatable.txt","1","reference.txt","1","2","y","temp");
    $parsing->vlookup("temp","1","index_file_ids.txt","2","3","y","final_downloadtable_$disease_abbr.txt");
    
    `mkdir $Analysispath/$disease_abbr/$tables` unless(-d "$Analysispath/$disease_abbr/$tables");
 
    if($Exp_Strategy eq "WGS")
    {
        `sort -k2,2 -k3,3 -k6,6 final_downloadtable_$disease_abbr.txt > final_downloadtable_$disease_abbr\_sorted.txt`;
        #No header in the output, thus no need to strip head!
        #pull_matched_tn_GDC(sorted bamlist,output file)
        $dwnld->pull_matched_tn_GDC("final_downloadtable_$disease_abbr\_sorted.txt","final_downloadtable_$disease_abbr\_$Exp_Strategy.txt");
        #filter table and remove bams aligning to NCBI36 or HG18;
        `cat final_downloadtable_$disease_abbr\_$Exp_Strategy.txt|grep NCBI36 -v | grep -v HG18 > WGS_tmp.txt;mv WGS_tmp.txt final_downloadtable_$disease_abbr\_$Exp_Strategy.txt`;
        copy("final_downloadtable_$disease_abbr\_$Exp_Strategy.txt","$Analysispath/$disease_abbr/$tables");
        #Downloads WGS BAMs and runs mpileups on them.
        #Dwld_WGSBam_and_do_mpileup(BamListfile,key directory,bam output directory,ref_fullpath,mpileup_outdir,user option(all,download or mpileups),line_num2split,"table directory")
        mkdir "$rna_wgs_dir/wgs_mpileups" unless(-d "$rna_wgs_dir/wgs_mpileups");
        $dwnld->Dwld_WGSBam_and_do_mpileup("$rna_wgs_dir/final_downloadtable_$disease_abbr\_$Exp_Strategy.txt","$rna_wgs_dir","$OUT_DIR","$database_path/GRCh37-lite.fa","$rna_wgs_dir/wgs_mpileups",$choice,$number,$VarScan_Path,$disease_abbr,"$database_path/hg19.fa","$Analysispath/$disease_abbr/$tables");
        copy("already_done_WGS.txt","$Analysispath/$disease_abbr/$tables")
    }
    elsif($Exp_Strategy eq "RNA-Seq")
    {
        #filter table and remove bams that may be aligned to HG18, HG19_Broad_variant, GRChr37-lite, NCBI36 or non aligned bams;
        `cat final_downloadtable_$disease_abbr.txt|grep NCBI36 -v|grep -v HG18 |grep -v HG19_Broad_variant |grep -v GRCh37-lite |grep -v NaN > final_downloadtable_$disease_abbr\_$Exp_Strategy.txt`;
  
	if ($choice ne "download")
	{
	     if(!(-d "$RNA_Path/$cds_sorted"))
	    {
		print STDERR "The directory $RNA_Path/$cds_sorted does not exist. It was moved, renamed or deleted\n";
		print STDERR "Please run script 1.3_Make_het_cds.pl.\n";
		exit;
	    }
	    
	    chdir $RNA_Path;
	    
	    `mkdir $ase` unless(-d "$ase");
	    
	    chdir "$ase";
	    
	    mkdir "$rna_mpileups" unless(-d "$rna_mpileups");
	    
	    chdir "$rna_wgs_dir";
	    
	    copy("final_downloadtable_$disease_abbr\_$Exp_Strategy.txt","$Analysispath/$disease_abbr/$tables");
	}
        else
	{
	    #Downloads RNA BAMs and runs mpileups on them.
	    #Dwld_RNASeq_Bam_and_do_mpileup(BamListfile,key directory,bam output directory,ref_fullpath,mpileup_outdir,path to cds_sorted directory,user option(all,download or mpileups),line_num2split)
        $dwnld->Dwld_RNASeq_Bam_and_do_mpileup("$rna_wgs_dir/$tables/final_downloadtable_$disease_abbr\_$Exp_Strategy.txt","$rna_wgs_dir","$OUT_DIR","$database_path/GRCh37-lite.fa","$ase/$rna_mpileups","$RNA_Path/$cds_sorted",$choice,$number);
	}
        
    }
}
else
{
    print "It seems that there is a table in the dir already:\n","$Analysispath/$disease_abbr/$tables/final_downloadtable_$disease_abbr\_$Exp_Strategy.txt\n";
    
    #Check gdc key and mv it to db first!
    #copyfile2newfullpath(path to gdc key,path where gdc key will be copied)
    $parsing->copyfile2newfullpath("$key","$rna_wgs_dir/gdc.key");
    
    if($Exp_Strategy eq "WGS")
    {
        #Dwld_WGSBam_and_do_mpileup(BamListfile,key directory,bam output directory,ref_fullpath,mpileup_outdir,user option(all,download or mpileups),line_num2split)
        mkdir "$rna_wgs_dir/wgs_mpileups" unless(-d "$rna_wgs_dir/wgs_mpileups");
        $dwnld->Dwld_WGSBam_and_do_mpileup("$Analysispath/$disease_abbr/$tables/final_downloadtable_$disease_abbr\_$Exp_Strategy.txt","$rna_wgs_dir","$OUT_DIR","$database_path/GRCh37-lite.fa","$rna_wgs_dir/wgs_mpileups",$choice,$number,$VarScan_Path,$disease_abbr,"$database_path/hg19.fa","$Analysispath/$disease_abbr/$tables");
        copy("already_done_WGS.txt","$Analysispath/$disease_abbr/$tables")
    }
    elsif($Exp_Strategy eq "RNA-Seq")
    {
	if ($choice ne "download")
	{
	    if(!(-d "$RNA_Path/$cds_sorted"))
	    {
		print STDERR "The directory $RNA_Path/$cds_sorted does not exist. It was moved, renamed or deleted\n";
		print STDERR "Please run script 1.3_Make_het_cds.pl.\n";
		exit;
	    }
	    
	    chdir $RNA_Path;

	    `mkdir $ase` unless(-d "$ase");
	    
	    mkdir "$ase/$rna_mpileups" unless(-d "$ase/$rna_mpileups");
	    
	    chdir "$rna_wgs_dir";
	}
	    #Dwld_RNASeq_Bam_and_do_mpileup(BamListfile,key directory,bam output directory,ref_fullpath,mpileup_outdir,path to cds_sorted directory,user option(all,download or mpileups),line_num2split)
	    $dwnld->Dwld_RNASeq_Bam_and_do_mpileup("$Analysispath/$disease_abbr/$tables/final_downloadtable_$disease_abbr\_$Exp_Strategy.txt","$rna_wgs_dir","$OUT_DIR","$database_path/GRCh37-lite.fa","$ase/$rna_mpileups","$RNA_Path/$cds_sorted",$choice,$number);   
    }
}

chdir "$rna_wgs_dir";

my @del_files = $parsing->get_only_files_in_dir("$rna_wgs_dir");
#Removes all files after the processes have finished.
for(my $i = 0;$i < scalar(@del_files);$i++)
{
    `rm "$del_files[$i]"`;
}

`rm -rf tmp*`;

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
