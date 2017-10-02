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

my $parsing = Driver_ASE_Lib::Parsing_Routines->new;
my $dwnld = Driver_ASE_Lib::Dwnld_WGS_RNA->new;

GetOptions(
    'cancer|c=s' => \my $cancer_type,#e.g. OV 
    'Expstrategy|E=s' => \my $Exp_Strategy,#e.g. Genotyping array
    'arraytype|a=s' =>\my $array_type,#e.g Genotypes
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("0_table");

if ($help)
{
    $parsing->usage("0_table");
}

if (!defined $cancer_type || !defined $Exp_Strategy || !defined $array_type)
{
    print STDERR "Cancer type, array type and/or experimental strategy was not entered!\n";
    $parsing->usage("0_table");
}

if ($Exp_Strategy eq "Genotyping array")
{
    if ("$array_type" ne "Genotypes" and "$array_type" ne "Copy number estimate")
    {
        print STDERR "Array type must be Genotypes or Copy number estimate, as those are the types that are used in this pipeline.\n";
        $parsing->usage("0_table");
    }
}
else
{
    print STDERR "The experimental strategy that was entered in was not the right one, it should be Genotyping array for this script.\n";
    $parsing->usage("0_table");
}

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my $Table_Dir = "$Analysispath/tables";
my $tables = "$cancer_type\_tables";
my $Geno_CNV_table = "$cancer_type.$array_type.id2uuid.txt";

$parsing->check_directory_existence("$database_path"); #check if directories or files exist

my @cancer;
if ($cancer_type =~ /,/)
{
    @cancer = split(",",$cancer_type);
}
else
{
   @cancer = $cancer_type; 
}

mkdir "$Driver_ASE_Dir/Analysis" unless(-d "$Driver_ASE_Dir/Analysis");

chdir $Analysispath;

`mkdir -p "$Table_Dir"`;

chdir $Table_Dir;

foreach my $cancer_type (@cancer)
{
    $cancer_type = $parsing->check_cancer_type($database_path,$cancer_type);
    next unless (defined $cancer_type);
    
    `mkdir $Analysispath/$cancer_type` unless(-d "$Analysispath/$cancer_type");
    
    if (!(-f "$Analysispath/$cancer_type/$tables/$Geno_CNV_table"))
    {
        #gets the manifest file from gdc and gets the UUIDs from it
        #gdc_parser(cancer type(e.g. OV),type of data (Genotyping array),data type)
        $dwnld->gdc_parser($cancer_type,$Exp_Strategy,$array_type);
        
        if ($array_type eq "Copy number estimate")
        {
            open (my $tangent,"$cancer_type.$array_type.result.txt");
            open (my $out,">$cancer_type\_tangent.txt");
            my @tanarray = <$tangent>;
            chomp(@tanarray);
            close ($tangent);
            @tanarray = grep{/tangent/}@tanarray;
            foreach (@tanarray)
            {
                print $out "$_\n";
            }
            close ($out);
            #metadata_collect(tangent.txt file,output file)
            $dwnld->metadata_collect("$cancer_type\_tangent.txt","$cancer_type\_$array_type\_Payload.txt");
        }
        #This filter only gets birdseed files. This is mainly for cancer types that have files that are not just birdseed
        elsif ($array_type eq "Genotypes")
        {
            open (my $birdseed,"$cancer_type.$array_type.result.txt");
            open (my $out,">$cancer_type\_birdseed.txt");
            my @birdseedarray = <$birdseed>;
            chomp(@birdseedarray);
            close ($birdseed);
            @birdseedarray = grep{/birdseed/}@birdseedarray;
            foreach (@birdseedarray)
            {
                print $out "$_\n";
            }
            close ($out);
            #puts the UUIDs in a payload file which will be used for the curl command
            #metadata_collect(birdseed.txt file,output file)
            $dwnld->metadata_collect("$cancer_type\_birdseed.txt","$cancer_type\_$array_type\_Payload.txt");
        }
        
        `curl --request POST --header \'Content-Type: application/json\' --data \@\'$cancer_type\_$array_type\_Payload.txt\' \'https://gdc-api.nci.nih.gov/legacy/files\' > \'$cancer_type\_$array_type.metadata.txt\'`;
        
        #matches UUID and TCGA ID
        #The columns of each cancer type may be different
        #QueryBase(.result.txt file from gdc_parser,query column,.metadata.txt file from the curl command,reqular expression,output file,column(s) to keep)
        $parsing->QueryBase("$cancer_type.$array_type.result.txt",1,"$cancer_type\_$array_type.metadata.txt",'TCGA-\w+-\w+-\w+',"t1.txt",0);
        $parsing->QueryBase("t1.txt",1,"$cancer_type\_$array_type.metadata.txt",'(tumor|blood|normal)',"$cancer_type.$array_type.id2uuid_query.txt",1);
	
	open (IDI,"$cancer_type.$array_type.id2uuid_query.txt");
	open (IDO,">$Geno_CNV_table");
	
	while (my $r = <IDI>)
	{
	    chomp($r);
	    my @a = split("\t",$r);
	    my $TCGA = $a[1];
	    my $sample = [split("-",$TCGA)]->[-1];
	    $TCGA =~ s/-[0-9]+[a-zA-Z]$//;
	    $sample =~ s/[a-zA-Z]$//;
            $sample =~s/^0//;
	    print IDO "$a[0]\t$a[1]\t$sample\t$a[2]\t$TCGA\n";
	}
	close (IDI);
	close (IDO);
        
        `mkdir $Analysispath/$cancer_type/$tables` unless(-d "$Analysispath/$cancer_type/$tables");
        
        copy("$Geno_CNV_table","$Analysispath/$cancer_type/$tables");
    }
    else
    {
        print "It seems that there is a table in the $Analysispath/$cancer_type/$tables directory already. Table: $Geno_CNV_table\n";
    }
}

chdir "$Analysispath";

`rm -rf "$Table_Dir"`;

print "All jobs have finished for $cancer_type.\n";
  
$time = localtime;
print "Script finished: $time.\n";

exit;
