#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$Bin/..";
use WGS_Analysis;
use Parsing_Routines;
use Cwd "realpath";
use File::Copy;
use Getopt::Long;
use strict;

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $wgs_analysis = Driver_ASE_Lib::WGS_Analysis->new;
my $parsing = Driver_ASE_Lib::Parsing_Routines->new;

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'readcutoff|r=i' => \my $read_cutoff,#e.g. 20
    'tfreq|t=i' => \my $tumor_freq,#e.g. 0.1
    'nfreq|n=i' => \my $normal_freq,
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage;

if($help)
{
    $parsing->usage("4.0");
}

if(!defined $disease_abbr)
{
    print "disease type was not entered!\n";
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

my $Analysispath = realpath("../../Analysis");
my $wgs_dwnlds = "wgs_dwnlds";
my $wgs_mpileups = "wgs_mpileups";
my $WGS_Path = "$Analysispath/$disease_abbr/WGS_Analysis";
my $ssoe_dir = "wgs_mpileups_ssoe_files";
my $finished_WGS = "$disease_abbr\_finished_analysis_WGS";
my $somatic = "somatic_variants";

#Checks if there is not Analysis directory.
if(!(-d "$Analysispath"))
{
    print STDERR "$Analysispath does not exist, it was either deleted, moved or renamed.\n";
    print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
    exit;
}
elsif(!(-d "$Analysispath/$disease_abbr"))
{
    print STDERR "$Analysispath/$disease_abbr does not exist, it was either deleted, moved or renamed.\n";
    print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
    exit;
}
elsif (!(-d "$Analysispath/$disease_abbr/$wgs_dwnlds"))
{
    print STDERR "$Analysispath/$disease_abbr/wgs_dwnlds does not exist. It was either deleted, moved or renamed.\n";
    print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
    exit;
}
elsif (!(-d "$Analysispath/$disease_abbr/$wgs_dwnlds/$wgs_mpileups"))
{
    print STDERR "$Analysispath/$disease_abbr/wgs_dwnlds/wgs_mpileups does not exist. It was either deleted, moved or renamed.\n";
    print STDERR "Please run script 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl.\n";
    exit;
}

`mkdir -p "$WGS_Path"` unless(-d "$WGS_Path");

my @ss = `ls $Analysispath/$disease_abbr/$wgs_dwnlds/$wgs_mpileups/`;
@ss = grep{/\.ss/}@ss;
my @o = `ls $Analysispath/$disease_abbr/$wgs_dwnlds/$wgs_mpileups/`;
@o = grep{/\.o/}@o;
my @e = `ls $Analysispath/$disease_abbr/$wgs_dwnlds/$wgs_mpileups/`;
@e = grep{/\.e/}@e;

if (@ss or @o or @e)
{
    `mkdir $WGS_Path/$ssoe_dir` unless(-d "$WGS_Path/$ssoe_dir");

    if(@ss)
    {
        foreach my $s(@ss)
        {
            chomp($s);
            `mv $Analysispath/$disease_abbr/$wgs_dwnlds/$wgs_mpileups/$s $WGS_Path/$ssoe_dir`;
        }
    }

    if(@o)
    {
        foreach my $out(@o)
        {
            chomp($out);
            `mv $Analysispath/$disease_abbr/$wgs_dwnlds/$wgs_mpileups/$out $WGS_Path/$ssoe_dir`;
        }   
    }
   
    if(@e)
    {
        foreach my $err(@e)
        {
            chomp($err);
            `mv $Analysispath/$disease_abbr/$wgs_dwnlds/$wgs_mpileups/$err $WGS_Path/$ssoe_dir`;   
        }
    }
}

mkdir "$Analysispath/$disease_abbr/$finished_WGS" unless (-d "$Analysispath/$disease_abbr/$finished_WGS");

chdir $WGS_Path;

if (-d "$WGS_Path/$ssoe_dir")
{
    if (-e "$ssoe_dir.tar.gz")
    {
        `tar -zcvkf $disease_abbr\_$ssoe_dir.tar.gz $ssoe_dir`;
        `mv $disease_abbr\_$ssoe_dir.tar.gz $Analysispath/$disease_abbr/$finished_WGS`;
    }
    else
    {
        `tar -zcvf $disease_abbr\_$ssoe_dir.tar.gz $ssoe_dir`;
        `mv $disease_abbr\_$ssoe_dir.tar.gz $Analysispath/$disease_abbr/$finished_WGS`;
    }
}

chdir "$Analysispath/$disease_abbr/$wgs_dwnlds";

print "Compressing wgs_mpileups directory\n";
`tar -zcvf $disease_abbr\_wgs_mpileups_varscan.tar.gz $wgs_mpileups`;
`mv $Analysispath/$disease_abbr/$wgs_dwnlds/$disease_abbr\_wgs_mpileups_varscan.tar.gz $Analysispath/$disease_abbr/$finished_WGS`;

chdir "$WGS_Path";

mkdir "$WGS_Path/$somatic";

`rm -f $WGS_Path/$somatic/*`;

#Varscan_filter(full path to wgs mpileups directory,full path to somatic_variants directory,readcutoff,VarType(e.g. Somatic),normal_alt_frq,tumor_alt_frq)
$wgs_analysis->Varscan_filter("$Analysispath/$disease_abbr/$wgs_dwnlds/$wgs_mpileups","$WGS_Path/$somatic",$read_cutoff,"Somatic",$normal_freq,$tumor_freq);

print "Compressing somatic_variants\n";
`tar -zcvf $disease_abbr\_$somatic.tar.gz $somatic`;
`mv $WGS_Path/$disease_abbr\_$somatic.tar.gz $Analysispath/$disease_abbr/$finished_WGS`;

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;
