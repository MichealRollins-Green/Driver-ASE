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

my $parsing = TCGA_Lib::Parsing_Routines->new;
my $dwnld = TCGA_Lib::Dwnld_WGS_RNA->new;

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV 
    'exp_strat|e=s' => \my $Exp_Strategy,#e.g. Genotyping array
    'array_type|a=s' =>\my $array_type,#e.g Genotypes
    'key|k=s'=> \my $key,
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("0");

if($help)
{
    $parsing->usage("0");
}

my $TCGA_Pipeline_Dir = realpath("../../");
mkdir "$TCGA_Pipeline_Dir/Analysis" unless(-d "$TCGA_Pipeline_Dir/Analysis");
my $Analysispath = realpath("../../Analysis");
my $Table_Dir = "$Analysispath/tables";
my $tables = "$disease_abbr\_tables";

if(!defined $disease_abbr || !defined $Exp_Strategy || !defined $array_type)
{
    print "disease type, array type and/or experimental strategy was not entered!\n";
    $parsing->usage("0");
}

if ($Exp_Strategy eq "Genotyping array")
{
    if("$array_type" ne "Genotypes" and "$array_type" ne "Copy number estimate")
    {
        print STDERR "data type must be Genotypes or Copy number estimate as those are the types that are used in this pipeline.\n";
        $parsing->usage("0");
    }
}
else
{
    print "The experimental strategy that was entered in was not the right one, it should be Genotyping array for this script.\n";
    $parsing->usage("0");
}

if(!defined $key or (!(-f $key)))
{
    print "gdc key fullpath was not entered or the fullpath to it was not correct!\n";
    print $key,"\n";
    $parsing->usage("0");
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

#Check if there is no Database Directory
if(!(-d "$TCGA_Pipeline_Dir/Database"))
{
    print STDERR "$TCGA_Pipeline_Dir/Database does not exist, it was either moved, renamed, deleted or has not been downloaded.\nPlease check the README.md file on the github page to find out where to get the Database directory.\n";
    exit;
}

chdir $Analysispath;

`mkdir -p "$Table_Dir"`;

chdir $Table_Dir or die "Can't change to directory $Table_Dir: $!\n";

foreach my $disease_abbr(@disease)
{
    `mkdir $Analysispath/$disease_abbr` unless(-d "$Analysispath/$disease_abbr");
    
    #Check gdc key and mv it to db first!
    #copyfile2newfullpath(path to gdc key file,path where gdc key file will be copied)
    $parsing->copyfile2newfullpath("$key","$Table_Dir/gdc.key");
    
    if(!(-f "$Analysispath/$disease_abbr/$tables/$disease_abbr.$array_type.id2uuid.txt"))
    {
        #gets the manifest file from gdc and gets the UUIDs from it
        #gdc_parser(cancer type(e.g. OV),type of data (Genotyping array),data type)
        $dwnld->gdc_parser($disease_abbr,$Exp_Strategy,$array_type);
        
        if($array_type eq "Copy number estimate")
        {
            open(my $tangent,"$disease_abbr.$array_type.result.txt") or die "Can't open file for input: $!\n";
            open(my $out,">$disease_abbr\_tangent.txt") or die "Can't open file for output: $!\n";
            my @tanarray = <$tangent>;
            chomp(@tanarray);
            close ($tangent);
            @tanarray = grep{/tangent/}@tanarray;
            foreach(@tanarray)
            {
                print $out "$_\n";
            }
            close($out);
            #metadata_collect(tangent.txt file,output file)
            $dwnld->metadata_collect("$disease_abbr\_tangent.txt","$disease_abbr\_$array_type\_Payload.txt");
        }
        #This filter only gets birdseed files. This is mainly for cancer types that have files that are not just birdseed
        elsif($array_type eq "Genotypes")
        {
            open(my $birdseed,"$disease_abbr.$array_type.result.txt") or die "Can't open $disease_abbr.$array_type.result.txt for input: $!\n";
            open(my $out,">$disease_abbr\_birdseed.txt") or die "Can't open file $disease_abbr\_birdseed.txt for output: $!\n";
            my @birdseedarray = <$birdseed>;
            chomp(@birdseedarray);
            close ($birdseed);
            @birdseedarray = grep{/birdseed/}@birdseedarray;
            foreach(@birdseedarray)
            {
                print $out "$_\n";
            }
            close($out);
            #puts the UUIDs in a payload file which will be used for the curl command
            #metadata_collect(birdseed.txt file,output file)
            $dwnld->metadata_collect("$disease_abbr\_birdseed.txt","$disease_abbr\_$array_type\_Payload.txt");
        }
        
        `curl --request POST --header \'Content-Type: application/json\' --data \@\'$disease_abbr\_$array_type\_Payload.txt\' \'https://gdc-api.nci.nih.gov/legacy/files\' > \'$disease_abbr\_$array_type.metadata.txt\'`;
        
        #matches UUID and TCGA ID
        #The columns of each cancer type may be different
        #QueryBase(.result.txt file from gdc_parser,query column,.metadata.txt file from the curl command,reqular expression,output file,column(s) to keep)
        $parsing->QueryBase("$disease_abbr.$array_type.result.txt",1,"$disease_abbr\_$array_type.metadata.txt",'TCGA-\w+-\w+-\w+',"t1.txt",0);
        $parsing->QueryBase("t1.txt",1,"$disease_abbr\_$array_type.metadata.txt",'(tumor|blood|normal)',"$disease_abbr.$array_type.id2uuid_query.txt",1);
	
	open(IDI,"$disease_abbr.$array_type.id2uuid_query.txt") or die "Can't open file $disease_abbr.$array_type.id2uuid_query.txt: $!\n";
	open(IDO,">$disease_abbr.$array_type.id2uuid.txt") or die "Can't open file $disease_abbr.$array_type.id2uuid.txt: $!\n";
	
	while(my $r = <IDI>)
	{
	    chomp($r);
	    my @a = split("\t",$r);
	    my $TCGA = substr($a[1],0,12);
	    my $sample = [split("-",$a[1])]->[-1];
	    $sample =~ s/\D+//;
            $sample =~s/^0//;
	    print IDO "$a[0]\t$a[1]\t$sample\t$a[2]\t$TCGA\n";
	}
	close(IDI);
	close(IDO);
        
        `mkdir $Analysispath/$disease_abbr/$tables` unless(-d "$Analysispath/$disease_abbr/$tables");
        
        copy("$disease_abbr.$array_type.id2uuid.txt","$Analysispath/$disease_abbr/$tables");
    }
    else
    {
        print "It seems that there is a table in the dir already: $disease_abbr.$array_type.id2uuid.txt\n";
    }
}

`rm -rf "$Table_Dir"`;

print "All jobs have finished for $disease_abbr.\n";
  
$time = localtime;
print "Script finished on $time.\n";

exit;
