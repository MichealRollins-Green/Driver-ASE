#!/usr/bin/perl -w

use MCE;
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

my $parsing = TCGA_Lib::Parsing_Routines->new;
my $dwnld = TCGA_Lib::Dwnld_WGS_RNA->new;

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV or OV,PRAD
    'exp_strat|e=s' => \my $Exp_Strategy,#e.g. WGS RNA-Seq
    'table_dir|b=s' => \my $Table_Dir,#directory where the table will be created.
    'key|k=s'=>\my $key,#path to the gdc.key.
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("3.0_table");
#If -h was specified on the command line then the usage of the program will be printed.
if($help)
{
    $parsing->usage("3.0_table");
}
#If the disease abbr or file type was not specified then an error will be printed and the usage of the program will be shown.
if(!defined $disease_abbr || !defined $Exp_Strategy)
{
    print "disease type and/or file type was not entered!\n";
    $parsing->usage("3.0_table");
}
#If the path to the gdc key was not specified then an error will be printed and the usage of the program will be shown.
if(!defined $key or (!( -f $key)))
{
    print "gdc key fullpath was not entered or the fullpath to it was not correct!\n";
    $parsing->usage("3.0_table");
}

if ($Exp_Strategy ne "RNA-Seq" and $Exp_Strategy ne "WGS")
{
    print STDERR "The file type entered must be RNA-Seq or WGS as these are what this pipeline deals with\n";
    exit;
}

my $TCGA_Pipeline_Dir = realpath("../../");

#Checks if there is not Database directory
if(!(-d "$TCGA_Pipeline_Dir/Database"))
{
    print STDERR "$TCGA_Pipeline_Dir/Database does not exist, it was either moved, deleted or has not been downloaded.\nPlease check the README.md file on the github page to find out where to get the Database directory.\n";
    exit;
}

my @disease;
if($disease_abbr =~ /,/)
{
    @disease = split(",",$disease_abbr);
}
else
{
   @disease = $disease_abbr; 
}

#Directory where all analysis data will be going in.
mkdir "$TCGA_Pipeline_Dir/Analysis" unless(-d "$TCGA_Pipeline_Dir/Analysis");
my $Analysispath=realpath("../../Analysis");

chdir "$Analysispath";

#defaults to a directory if no output directory was specified in the command line.
if(!defined $Table_Dir)
{
    if($Exp_Strategy eq "RNA-Seq")
    {
        $Table_Dir = "$Analysispath/tables";
    }
    elsif($Exp_Strategy eq "WGS")
    {
        $Table_Dir = "$Analysispath/tables";
    }
    else
    {
        print STDERR "File type must be either RNA-Seq or WGS.\n";
        exit;
    }
    
    print STDERR "no output directory entered, using $Table_Dir as default\n";
    `mkdir -p $Table_Dir` unless(-d $Table_Dir);
    }
 else
 {
    unless(-d $Table_Dir)
    {      
        print STDERR "Please provide fullpath for output dir: \nThe dir does not exist: $Table_Dir\n";
        exit;
    }
    $Table_Dir =~ s/\/$//;
    if($Exp_Strategy eq "RNA-Seq")
    {   
        `mkdir -p "$Table_Dir/tables"` unless(-d "$Table_Dir/tables");
        $Table_Dir .= "/tables";
    }
    elsif($Exp_Strategy eq "WGS")
    {
        `mkdir -p "$Table_Dir/tables"` unless(-d "$Table_Dir/tables");
        $Table_Dir .= "/tables";
    }
 } 
    
chdir $Table_Dir;

foreach my $disease_abbr(@disease)
{
    #Makes the directory where all of the results will be processed and stored.
    `mkdir $Analysispath/$disease_abbr` unless(-d "$Analysispath/$disease_abbr");
    #Check gdc key and mv it to db first!
    #copyfile2newfullpath(path to gdc key,path where gdc key will be copied)
    $parsing->copyfile2newfullpath("$key","$Table_Dir/gdc.key");
    
    #Checks if a table file does not exist in the tables directory of the cancer type.
    if(!(-f "$Analysispath/$disease_abbr/$disease_abbr\_tables/final_downloadtable_$disease_abbr\_$Exp_Strategy.txt"))
    {
        #parses the gdc website manifest of the specified cancer type and prints to a results file.
        #gdc_parser(cancer type(e.g. OV),type(RNA-Seq or WGS))
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
        
        `curl --request POST --header \'Content-Type: application/json\' --data \@Payload.txt 'https://gdc-api.nci.nih.gov/legacy/files' > $disease_abbr\_metadata.txt`;
        
        open(my $meta,"$disease_abbr\_metadata.txt") or die "Can\'t open file for input: $!";
        open my $ME,">$disease_abbr.edit.metadata.txt" or die "Can\'t open file for output: $!";
        chomp(my @metaedit = <$meta>);
        close($meta);
        for(my $i = 0;$i < scalar(@metaedit);$i++)
        {
            $metaedit[$i]=~s/\r//g;
        }
        @metaedit = grep{!/\t\s+\t/}@metaedit;
        print $ME "$_\n" for @metaedit;
        close($ME);
        
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
        $dwnld->ref_parse("meta_ids.txt","$Table_Dir","reference.txt");
        
        #index_ids(.result.txt,output file)
        $dwnld->index_ids("$disease_abbr.result.txt","Payload.txt");
        
        `curl --request POST --header \'Content-Type: application/json\' --data \@Payload.txt 'https://gdc-api.nci.nih.gov/legacy/files' > index_file_ids.txt`;
        
       #vlookup(lookupFile,queryCol,sourceFile,lookupCol,returnCol(s),append(y/n),outputFile)
    
        #e.g. vlookup(lookupfile,3,sourcefile,4,"1,2,4,6","y",outputFile)

        #Will search each column 3 entry of lookupfile within column 4 of sourceFile and will append columns 1,2,4,6 of sourceFile to the end of each row in lookupfile.
        
        #N.B. only works on tab-delimited files
        $parsing->vlookup("$disease_abbr.datatable.txt","1","reference.txt","1","2","y","temp");
        $parsing->vlookup("temp","1","index_file_ids.txt","2","3","y","final_downloadtable_$disease_abbr.txt");
        
        `mkdir $Analysispath/$disease_abbr/$disease_abbr\_tables` unless(-d "$Analysispath/$disease_abbr/$disease_abbr\_tables");
        
        if($Exp_Strategy eq 'WGS')
        {
            `sort -k2,2 -k3,3 -k6,6 final_downloadtable_$disease_abbr.txt > final_downloadtable_$disease_abbr\_sorted.txt`;
            #No header in the output, thus no need to strip head!
            #pull_matched_tn_GDC(sorted bamlist,output file)
            $dwnld->pull_matched_tn_GDC("final_downloadtable_$disease_abbr.txt","final_downloadtable_$disease_abbr\_$Exp_Strategy.txt");
            #filter table and remove bams aligning to NCBI36 or HG18;
            `cat final_downloadtable_$disease_abbr\_$Exp_Strategy.txt|grep NCBI36 -v | grep -v HG18 > WGS_tmp.txt;mv WGS_tmp.txt final_downloadtable_$disease_abbr\_$Exp_Strategy.txt`;
            copy("final_downloadtable_$disease_abbr\_$Exp_Strategy.txt","$Analysispath/$disease_abbr/$disease_abbr\_tables");
        }
        elsif($Exp_Strategy eq "RNA-Seq")
        {
           `cat final_downloadtable_$disease_abbr.txt |grep NCBI36 -v|grep -v HG18 |grep -v HG19_Broad_variant |grep -v GRCh37-lite |grep -v NaN > final_downloadtable_$disease_abbr\_sort.txt`;
           `sort -k2,2 final_downloadtable_$disease_abbr\_sort.txt > final_downloadtable_$disease_abbr\_$Exp_Strategy.txt`;     
            copy("final_downloadtable_$disease_abbr\_$Exp_Strategy.txt","$Analysispath/$disease_abbr/$disease_abbr\_tables");
        }
    }
    else
    {
       print "It seems that there is a table in the dir already: final_downloadtable_$disease_abbr\_$Exp_Strategy.txt\n";
    }
}

`rm -rf $Table_Dir`;

print "All jobs have finished for $disease_abbr.\n";
  
$time = localtime;
print "Script finished on $time.\n";

exit;
