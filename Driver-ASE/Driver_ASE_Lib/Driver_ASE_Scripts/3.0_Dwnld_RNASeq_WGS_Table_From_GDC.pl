#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$Bin/..";
use Parsing_Routines;
use Dwnld_WGS_RNA;
use File::Copy;
use Getopt::Long;
use Cwd 'realpath';
use File::Basename;
no warnings 'once';
use autodie;

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $parsing = Driver_ASE_Lib::Parsing_Routines->new;
my $dwnld = Driver_ASE_Lib::Dwnld_WGS_RNA->new;

GetOptions(
    'cancer|c=s' => \my $cancer_type,#e.g. OV or OV,PRAD
    'Expstrategy|E=s' => \my $Exp_Strategy,#e.g. WGS RNA-Seq
    'download|d=s' => \my $dwnld_cmd,#curl or aria2c (if aria or aria2 is entered, it changes them to aria2c as that is the command)
    'overlap|o=s' => \my $overlap, #choice to overlap IDs
    'key|k=s' => \my $key,#path to gdc key
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("3.0_table");

#If -h was specified on the command line then the usage of the program will be printed.
if ($help)
{
    $parsing->usage("3.0_table");
}

#If the disease abbr or file type was not specified then an error will be printed and the usage of the program will be shown.
if (!defined $cancer_type || !defined $Exp_Strategy)
{
    print STDERR "Cancer type and/or experimental strategy was not entered!\n";
    $parsing->usage("3.0_table");
}

if ($Exp_Strategy ne "RNA-Seq" and $Exp_Strategy ne "WGS")
{
    print STDERR "The experimental strategy entered must be RNA-Seq or WGS as these are what this pipeline deals with!\n";
    exit;
}

my $Driver_ASE_Dir = realpath("../../");
my $Table_Dir = "tables";
my $tables = "$cancer_type\_tables";
my $database_path = "$Driver_ASE_Dir/Database";
my $Intersect = "$cancer_type\_RNA_WGS_Geno.txt";
my $downloadtable = "final_downloadtable_$cancer_type";
my $RNA_table = "$downloadtable\_RNA-Seq.txt";
my $WGS_table = "$downloadtable\_WGS.txt";
my $Genotypes_table = "$cancer_type.Genotypes.id2uuid.txt";
my $RNA_table_overlap = "$downloadtable\_RNA-Seq_overlap.txt";
my $WGS_table_overlap = "$downloadtable\_WGS_overlap.txt";
my $Genotypes_table_overlap = "$cancer_type.Genotypes.id2uuid_overlap.txt";
my $reference_file = "$cancer_type\_$Exp_Strategy\_reference.txt";

#If the path to the gdc key was not specified then an error will be printed and the usage of the program will be shown.
if (!defined $key)
{
    print STDERR "GDC key fullpath was not entered or the fullpath to it was not correct!\n";
    $parsing->usage("3.0_table");
}

$parsing->check_directory_existence("$database_path","$key"); #check if directories or files exist

#Defaults to curl if no download command was specified
if (!defined $dwnld_cmd)
{
    $dwnld_cmd = "curl";
    print "No download command specified! defaulting to $dwnld_cmd\n";
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
    exit;
}

$parsing->CheckSoftware("$dwnld_cmd");

#Directory where all analysis data will be going in.
mkdir "$Driver_ASE_Dir/Analysis" unless(-d "$Driver_ASE_Dir/Analysis");
my $Analysispath = realpath("../../Analysis");

if (!defined $overlap)
{
    $overlap = "no";
}
elsif (lc $overlap ne "n" and lc $overlap ne "no" and lc $overlap ne "y" and lc $overlap ne "yes")
{
    print STDERR "Incorrect option entered! options are n|no|y|yes.\n";
    print $parsing->usage("3_table");
}

chdir "$Analysispath";

#defaults to a directory if no output directory was specified in the command line.
my $Input_Dir;
if ($Exp_Strategy eq "RNA-Seq")
{   
   `mkdir -p "$Table_Dir/rna"` unless(-d "$Table_Dir/rna");
   $Input_Dir = "$Table_Dir";
   $Table_Dir .= "/rna";
}
elsif ($Exp_Strategy eq "WGS")
{
   `mkdir -p "$Table_Dir/wgs"` unless(-d "$Table_Dir/wgs");
   $Input_Dir = "$Table_Dir";
   $Table_Dir .= "/wgs";
}
else
{
    print STDERR "Experimental strategy must be either RNA-Seq or WGS!\n";
    exit;
}

$Table_Dir = realpath("$Table_Dir");

my @cancer;
if ($cancer_type =~ /,/)
{
    @cancer = split(",",$cancer_type);
}
else
{
   @cancer = $cancer_type; 
}

foreach my $cancer_type (@cancer)
{
    $cancer_type = $parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid
    next unless (defined $cancer_type);

    #Makes the directory where all of the results will be processed and stored.
    `mkdir "$Analysispath/$cancer_type"` unless(-d "$Analysispath/$cancer_type");
 
    chdir "$Table_Dir";
    
    #Check gdc key and mv it to db first!
    #copyfile2newfullpath(path to gdc key file,path where gdc key file will be copied)
    $parsing->copyfile2newfullpath("$key","$Table_Dir/gdc.key");
    
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
        `curl --request POST --header "Content-Type: application/json" --data \@Payload.txt 'https://api.gdc.cancer.gov/v0/legacy/files' > $cancer_type\_metadata.txt`;
        
        #Uses the metadata file to get data(i.e. UUID, TCGA ID) an prints the data for the UUIDs to a datatable file.
        #parse_patient_id(_metadata.txt file created from curl,output file)
        $dwnld->parse_patient_id("$cancer_type\_metadata.txt","$cancer_type.datatable.txt");
        
        #metadata_ids(.result.txt,output file)
        $dwnld->metadata_ids("$cancer_type.result.txt","Payload.txt");
        
        `curl --request POST --header \'Content-Type: application/json\' --data \@Payload.txt 'https://api.gdc.cancer.gov/v0/legacy/files' > $cancer_type\_metadata.txt`;
     
        open (my $meta,"$cancer_type\_metadata.txt");
        open (my $ME,">$cancer_type.edit.metadata.txt");
        chomp(my @metaedit = <$meta>);
        close ($meta);
        for (my $i = 0;$i < scalar(@metaedit);$i++)
        {
        $metaedit[$i] =~ s/\r//g;
        }
        @metaedit = grep{!/\t\s+\t/}@metaedit;
        print $ME "$_\n" for @metaedit;
        close ($ME);
        
        #pulls the UUIDs from the edit metadata file
        #pull_column(.edit.metadata.txt,column(s) to pull,output file)
        $parsing->pull_column("$cancer_type.edit.metadata.txt","8","temp.UUID");
        
        #Strips the headers in the file.
        #strip_head(file to strip headers,output file)
        $parsing->strip_head("temp.UUID","$cancer_type.UUID.txt");
        `rm temp.UUID`;
        
        #parse_meta_id(.edit.metadata.txt,output file from strip_head,output file)
        $dwnld->parse_meta_id("$cancer_type.edit.metadata.txt","$cancer_type.UUID.txt","meta_ids.txt");
        
        #ref_parse(output file from parse_meta_id,directory building table,output file)
        $dwnld->ref_parse("meta_ids.txt","$Table_Dir","$reference_file",$dwnld_cmd);
        
        #index_ids(.result.txt,output file)
        $dwnld->index_ids("$cancer_type.result.txt","Payload.txt");
        
        `curl --request POST --header \'Content-Type: application/json\' --data \@Payload.txt 'https://api.gdc.cancer.gov/v0/legacy/files' > index_file_ids.txt`;
       
        #vlookup(lookupFile,queryCol,sourceFile,lookupCol,returnCol(s),append(y/n),outputFile)
        #e.g. vlookup(lookupfile,3,sourcefile,4,"1,2,4,6","y",outputFile)
        #Will search each column 3 entry of lookupfile within colum 4 of sourceFile and will append columns 1,2,4,6 of sourceFile to the end of each row in lookupfile.
        #N.B. only works on tab-delimited files
        $parsing->vlookup("$cancer_type.datatable.txt","1","$reference_file","1","2","y","temp");
        
        if ($Exp_Strategy eq "RNA-Seq")
        {
            $parsing->vlookup("temp","1","index_file_ids.txt","3","1","y","$downloadtable.txt");
        }
        else
        {
            $parsing->vlookup("temp","1","index_file_ids.txt","3","2","y","$downloadtable.txt");
        }
        
        `mkdir "$Analysispath/$cancer_type/$tables"` unless(-d "$Analysispath/$cancer_type/$tables");
        
        if ($Exp_Strategy eq 'WGS')
        {
            `sort -k2,2 -k3,3 -k6,6 $downloadtable.txt > $downloadtable\_sorted.txt`;
            #No header in the output, thus no need to strip head!
            #pull_matched_tn_GDC($downloadtable_sorted.txt file,output file)
            $dwnld->pull_matched_tn_GDC("$downloadtable\_sorted.txt","$WGS_table");
            $parsing->vlookup("$WGS_table",1,"$cancer_type.result.txt",1,4,"y","$downloadtable\_$Exp_Strategy\_size.txt");
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
                print SO $size_wgs,"\n";
            }
            close (SIZE);
            close (SO);
            
            #filter table and remove bams aligning to NCBI36 or HG18;
            `cat $downloadtable\_$Exp_Strategy\_convert.txt | grep -v NCBI36| grep -v HG18 > WGS_tmp.txt;mv WGS_tmp.txt $WGS_table`;
           copy("$WGS_table","$Analysispath/$cancer_type/$tables");
        }
        elsif ($Exp_Strategy eq "RNA-Seq")
        {
           `cat $downloadtable.txt |grep -v NCBI36|grep -v HG18 |grep -v HG18_Broad_variant|grep -v HG19_Broad_variant |grep -v GRCh37-lite > $downloadtable\_sort.txt`;
           `sort -k2,2 $downloadtable\_sort.txt > $RNA_table`;     
            copy("$RNA_table","$Analysispath/$cancer_type/$tables");
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
}

if (lc $overlap eq "y" || lc $overlap eq "yes")
{
    $parsing->check_directory_existence("$Analysispath/$cancer_type/$tables/$RNA_table","$Analysispath/$cancer_type/$tables/$WGS_table","$Analysispath/$cancer_type/$tables/$Genotypes_table"); #checks if RNA-Seq, WGS and Genotypes tables exist before overlapping them
    if (!(-s "$Analysispath/$cancer_type/$tables/$RNA_table" == 0) and !(-s "$Analysispath/$cancer_type/$tables/$WGS_table" == 0) and !(-s "$Analysispath/$cancer_type/$tables/$Genotypes_table" == 0))
    {
        #Overlap_RNA_WGS_Geno($Analysispath/$cancer_type/$tables (path to cancer type table directory),$RNA_table (RNA-Seq table file),$WGS_table (WGS_table file),$Genotypes_table (Genotypes table file),$Intersect (file to output overlapped results),$Exp_Strategy (RNA-Seq or WGS),$cancer_type (e.g. OV))
        $parsing->Overlap_RNA_WGS_Geno("$Analysispath/$cancer_type/$tables","$RNA_table","$WGS_table","$Genotypes_table","$RNA_table_overlap","$WGS_table_overlap","$Genotypes_table_overlap","$Intersect","$Exp_Strategy","$cancer_type");
    }
    else
    {
        print "Either $RNA_table, $WGS_table and/or $Genotypes_table is empty!\n";
        exit;
    }
}

chdir $Analysispath;
`rm -rf "$Input_Dir"`;

print "All jobs have finished for $cancer_type.\n";
  
$time = localtime;
print "Script finished on $time.\n";

exit;
