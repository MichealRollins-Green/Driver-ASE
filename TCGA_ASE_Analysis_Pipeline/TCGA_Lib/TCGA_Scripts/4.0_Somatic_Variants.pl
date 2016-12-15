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

my $wgs_analysis = TCGA_Lib::WGS_Analysis->new;
my $parsing = TCGA_Lib::Parsing_Routines->new;

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'wgs_dir|w=s' => \my $wgs_dir,#directory that holds the wgs_mpileups directory.
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("4.0");

if($help)
{
    $parsing->usage("4.0");
}

if(!defined $disease_abbr || !defined $wgs_dir)
{
    print "disease type and/or wgs directory was not entered!\n";
    $parsing->usage("4.0");
}

my $Analysispath = realpath("../../Analysis");

#Checks if there is not Analysis directory.
if(!(-d "$Analysispath"))
{
    print STDERR "$Analysispath does not exist, it was either deleted, moved or the script that creates it wasn't ran.\n";
    exit;
}
elsif(!(-d "$Analysispath/$disease_abbr"))
{
    print STDERR "$Analysispath/$disease_abbr does not exist, it was either deleted, moved or the script that creates it wasn't ran.\n";
    exit;
}
elsif (!(-d "$Analysispath/$disease_abbr/$wgs_dir"))
{
    print STDERR "$Analysispath/$disease_abbr/$wgs_dir does not exist. It was either deleted, moved or no bams have been processed with mpileups and varscan.\n";
    exit;
}
elsif (!(-d "$Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups"))
{
    print STDERR "$Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups does not exist. It was either deleted, moved or no bams have been processed with mpileups and varscan.\n";
    exit;
}

my $WGS_Path = "$Analysispath/$disease_abbr/WGS_Analysis";

`mkdir -p "$WGS_Path"` unless(-d "$WGS_Path");

my @ss = `ls $Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups/`;
@ss = grep{/\.ss/}@ss;
my @o = `ls $Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups/`;
@o = grep{/\.o/}@o;
my @e = `ls $Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups/`;
@e = grep{/\.e/}@e;

if(@ss)
{
    `mkdir $WGS_Path/wgs_mpileups_ssoe_files` unless(-d "$WGS_Path/wgs_mpileups_ssoe_files");
    foreach my $s(@ss)
    {
        chomp($s);
        `mv $Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups/$s $WGS_Path/wgs_mpileups_ssoe_files`;
    }
}

if(@o)
{
    `mkdir $WGS_Path/wgs_mpileups_ssoe_files` unless(-d "$WGS_Path/wgs_mpileups_ssoe_files");
    foreach my $out(@o)
    {
        chomp($out);
        `mv $Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups/$out $WGS_Path/wgs_mpileups_ssoe_files`;
    }   
}
   
if(@e)
{
    `mkdir $WGS_Path/wgs_mpileups_ssoe_files` unless(-d "$WGS_Path/wgs_mpileups_ssoe_files");
    foreach my $err(@e)
    {
        chomp($err);
        `mv $Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups/$err $WGS_Path/wgs_mpileups_ssoe_files`;   
    }
}

mkdir "$Analysispath/$disease_abbr/$disease_abbr\_finished_analysis_WGS";

chdir $WGS_Path;

if (-d "$WGS_Path/wgs_mpileups_ssoe_files")
{
    if (-e "wgs_mpileups_ssoe_files.tar.gz")
    {
        `tar -zcvkf wgs_mpileups_ssoe_files.tar.gz wgs_mpileups_ssoe_files`;
        `mv wgs_mpileups_ssoe_files.tar.gz $Analysispath/$disease_abbr/$disease_abbr\_finished_analysis_WGS`;
    }
    else
    {
        `tar -zcvf wgs_mpileups_ssoe_files.tar.gz wgs_mpileups_ssoe_files`;
        `mv wgs_mpileups_ssoe_files.tar.gz $Analysispath/$disease_abbr/$disease_abbr\_finished_analysis_WGS`;
    }
}

chdir "$Analysispath/$disease_abbr/$wgs_dir";

print "Compressing wgs_mpileups directory\n";
`tar -zcvf wgs_mpileups_varscan.tar.gz wgs_mpileups`;
mkdir "$Analysispath/$disease_abbr/$disease_abbr\_finished_analysis_WGS" unless(-d "$Analysispath/$disease_abbr/$disease_abbr\_finished_analysis_WGS");
`mv $Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups_varscan.tar.gz $Analysispath/$disease_abbr/$disease_abbr\_finished_analysis_WGS`;

chdir "$WGS_Path";

mkdir "$WGS_Path/somatic_variants";

`rm -f $WGS_Path/somatic_variants/*`;

#Varscan_filter(full path to wgs mpileups directory,full path to somatic_variants directory,readcutoff,VarType(e.g. Somatic),normal_alt_frq,tumor_alt_frq)
$wgs_analysis->Varscan_filter("$Analysispath/$disease_abbr/$wgs_dir/wgs_mpileups","$WGS_Path/somatic_variants",20,"Somatic",0.1,0.1);

print "Compressing somatic_variants\n";
`tar -zcvf somatic_variants.tar.gz somatic_variants`;
`mv $WGS_Path/somatic_variants.tar.gz $Analysispath/$disease_abbr/$disease_abbr\_finished_analysis_WGS`;

print "All jobs have finished for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;