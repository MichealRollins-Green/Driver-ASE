#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$Bin/..";
use WGS_Analysis;
use Parsing_Routines;
use Cwd "realpath";
use Getopt::Long;
use strict;
use autodie;
no warnings 'once';

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $wgs_analysis = Driver_ASE_Lib::WGS_Analysis->new;
my $parsing = Driver_ASE_Lib::Parsing_Routines->new;

GetOptions(
    'cancer|c=s' => \my $cancer_type, #e.g. OV
    'readcutoff|r=i' => \my $read_cutoff, #e.g. 20
    'tfreq|t=i' => \my $tumor_freq, #e.g. 0.1
    'nfreq|n=i' => \my $normal_freq,
    'overlap|o=s' => \my $overlap,
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("4.0");

if ($help)
{
    $parsing->usage("4.0");
}

if (!defined $cancer_type)
{
    print STDERR "Cancer type was not entered!\n";
    $parsing->usage("4.0");
}

#default to no if option to overlap ids was not specified in the command line
if (!defined $overlap)
{
    $overlap = "no";
}
elsif (lc $overlap ne "y" and lc $overlap ne "yes" and lc $overlap ne "n" and lc $overlap ne "no")
{
    print STDERR "Incorrect option specified! Must be yes|y|no|n.\n";
    $parsing->usage("4.0");
}

if (!defined $read_cutoff)
{
    $read_cutoff = 20;
}

if (!defined $tumor_freq)
{
    $tumor_freq = 0.1;
}

if (!defined $normal_freq)
{
    $normal_freq = 0.1;
}

my $Driver_ASE_Dir = realpath("../../");
my $database_path = "$Driver_ASE_Dir/Database";
my $Analysispath = realpath("../../Analysis");
my ($wgs_dwnlds,$wgs_mpileups,$somatic) = ("wgs_dwnlds","wgs_mpileups","somatic_variants");
my $WGS_Path = "$Analysispath/$cancer_type/WGS_Analysis";
my ($RNA_table_file,$RNA_table_overlap) = ("final_downloadtable_$cancer_type\_RNA-Seq.txt","final_downloadtable_$cancer_type\_RNA-Seq_overlap.txt");
my ($WGS_table_file,$WGS_table_overlap) = ("final_downloadtable_$cancer_type\_WGS.txt","final_downloadtable_$cancer_type\_WGS_overlap.txt");
my ($Geno_table_file,$Geno_table_overlap) = ("$cancer_type.Genotypes.id2uuid.txt","$cancer_type.Genotypes.id2uuid_overlap.txt");
my $Intersect = "$cancer_type\_RNA_WGS_Geno.txt";
my $tables = "$cancer_type\_tables";
my $Exp_Strategy = "WGS";
my ($WGS_Mpileups,$WGS_Mpileups_Overlapped) = ("WGS_Mpileups_Full.txt","WGS_Mpileups_Overlap.txt");
my $somatic_not_done = "$WGS_Path/not_done_somatic.txt";

$parsing->check_directory_existence("$database_path","$Analysispath","$Analysispath/$cancer_type","$Analysispath/$cancer_type/$wgs_dwnlds","$Analysispath/$cancer_type/$wgs_dwnlds/$wgs_mpileups"); #check if directories or files exist

$parsing->check_cancer_type($database_path,$cancer_type); #checks if the cancer type entered is valid

`mkdir -p "$WGS_Path"` unless(-d "$WGS_Path");

chdir "$WGS_Path";

mkdir "$WGS_Path/$somatic" unless (-d "$WGS_Path/$somatic");

#check if user wants to overlap RNA-Seq,WGS and Genotypes IDs
if (lc $overlap eq "y" || lc $overlap eq "yes")
{
    $parsing->check_directory_existence("$Analysispath/$cancer_type/$tables/$Geno_table_file","$Analysispath/$cancer_type/$tables/$RNA_table_file","$Analysispath/$cancer_type/$tables/$WGS_table_file");
    
    if (!(-s "$Analysispath/$cancer_type/$tables/$RNA_table_file" == 0) and !(-s "$Analysispath/$cancer_type/$tables/$WGS_table_file" == 0) and !(-s "$Analysispath/$cancer_type/$tables/$Geno_table_file" == 0))
    {
        #Overlap_RNA_WGS_Geno($Analysispath/$cancer_type/$tables (path to directory that holds the table files for RNA-Seq, WGS and Genotypes), $RNA_table_file (RNA-Seq table file), $WGS_table_file (WGS table file), $Geno_table_file (Genotypes table file), $Intersect (file that holds intersected data from RNA-Seq, WGS and Genotypes tables), $Exp_Strategy (Experimental Strategy (e.g. WGS)), $cancer_type (the cancer type (e.g. PRAD)))
        $parsing->Overlap_RNA_WGS_Geno("$Analysispath/$cancer_type/$tables","$RNA_table_file","$WGS_table_file","$Geno_table_file","$RNA_table_overlap","$WGS_table_overlap","$Geno_table_overlap","$Intersect","$Exp_Strategy","$cancer_type");
        chdir "$WGS_Path"; #changes back to this directory after overlapping as the Overlap_RNA_WGS_Geno routine changes to the $Analysispath/$cancer_type/$tables directory
        `ls $Analysispath/$cancer_type/$wgs_dwnlds/$wgs_mpileups > wgs_mpileups.txt`;
        
        open (MPO,">wgs_look.txt");
        my @mp = `cat wgs_mpileups.txt`;
        foreach my $mps (@mp)
        {
            chomp($mps);
            my @mpsplit = split("\\.",$mps); #split file names by . (dot)
            print MPO "$mpsplit[0]\t$mpsplit[1]\t$mpsplit[2]\t$mpsplit[3]\n";
        }
        close (MPO);
        $parsing->vlookup("wgs_look.txt",2,"$Analysispath/$cancer_type/$tables/$WGS_table_overlap",2,2,"y","wgs_mpileups_NaN.txt"); #new $WGS_table file contains only overlap IDs and matches those IDs to the IDs in the wgs_look.txt file
        my @filt_mpileups = `cat wgs_mpileups_NaN.txt`;
        open (MPO,">$WGS_Mpileups_Overlapped");
        foreach my $mpm (@filt_mpileups)
        {
            chomp($mpm);
            if ($mpm !~ /\tNaN+$/i)
            {
                my @new_mp = split("\t",$mpm);
                pop @new_mp; #pops the last element as it is just the IDs that were matched and printed from the vlookup routine
                my $joined_mp = join(".",$new_mp[0],$new_mp[1]); #re-joins the file name after being split and used in vlookup
                my $joined_mpm = join(".",@new_mp);
                print MPO "$Analysispath/$cancer_type/$wgs_dwnlds/$wgs_mpileups/$joined_mpm\t$joined_mp\t$Analysispath/$cancer_type/$wgs_dwnlds\t$wgs_mpileups/$joined_mpm\n"; #prints the path to the mpileup file
            }
        }
        close (MPO);      
    }
    else
    {
        print "Either $RNA_table_file, $WGS_table_file and/or $Geno_table_file is empty!\n";
        exit;
    }
}

`ls $WGS_Path/$somatic > $WGS_Path/$somatic\_done.txt`;

print "Running Varscan_filter.\n";
#Varscan_filter($Analysispath/$cancer_type/$wgs_dwnlds/$wgs_mpileups (path to ), $WGS_Path/$somatic (directory where somatic variant data will be stored), $read_cutoff (readcutoff e.g. 20),VarType(e.g. Somatic), $normal_freq (normal alternate frequency e.g. 0.1), $tumor_freq (tumor alternate frequency e.g. 0.1))
if (lc $overlap eq "y" or lc $overlap eq "yes")
{
    $wgs_analysis->filter_not_done_somatic("$somatic_not_done","$WGS_Path/$WGS_Mpileups_Overlapped","$WGS_Path/$somatic\_done.txt","$WGS_Path","$somatic");
    $wgs_analysis->Varscan_filter("$somatic_not_done","$WGS_Path/$somatic",$read_cutoff,"Somatic",$normal_freq,$tumor_freq);
}
else
{
    my @mpile = `ls $Analysispath/$cancer_type/$wgs_dwnlds/$wgs_mpileups`;
    open (MPO,">$WGS_Mpileups");
    
    for my $mp (@mpile)
    {
        chomp($mp);
        my @sp_mp = split("\\.",$mp);
        my $join_mp = join(".",$sp_mp[0],$sp_mp[1]); #join UUID and TCGA ID
        print MPO "$Analysispath/$cancer_type/$wgs_dwnlds/$wgs_mpileups/$mp\t$join_mp\t$Analysispath/$cancer_type/$wgs_dwnlds\t$wgs_mpileups/$mp\n";
    }
    close (MPO);
    
    $wgs_analysis->filter_not_done_somatic("$somatic_not_done","$WGS_Path/$WGS_Mpileups","$WGS_Path/$somatic\_done.txt","$WGS_Path","$somatic");
    $wgs_analysis->Varscan_filter("$somatic_not_done","$WGS_Path/$somatic",$read_cutoff,"Somatic",$normal_freq,$tumor_freq);
}

print "All jobs have finished for $cancer_type.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
