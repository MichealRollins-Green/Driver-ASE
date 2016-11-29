#!/usr/bin/perl -w
BEGIN{
	    push @INC, "/home/hpc3256/NGS_lib/Linux_codes_SAM/usr/share/perl5";
	}
use MCE::Map;
use FindBin qw($Bin);
use lib "$Bin/..";
use Imputation_Plink;
use Parsing_Routines;
use Cwd 'realpath';
use Cwd;
use Getopt::Long;
use strict;

my $time = localtime;
print "Script started on $time.\n";

#Changes to the directory of the script executing;
chdir $Bin;

my $impute_plink = TCGA_Lib::Imputation_Plink->new;
my $parsing = TCGA_Lib::Parsing_Routines->new;
my $TCGA_Pipeline_Dir = realpath("../../");
my $Analysispath = realpath("../../Analysis");

GetOptions(
    'disease|d=s' => \my $disease_abbr,#e.g. OV
    'output|o=s' => \my $imputation,
    'help|h' => \my $help
) or die "Incorrect options!\n",$parsing->usage("1.2");

if($help)
{
    $parsing->usage("1.2");
}

if(!defined $disease_abbr)
{
    print "disease type was not entered!\n";
    $parsing->usage("1.2");
}

my $RNA_Path= "$Analysispath/$disease_abbr/RNA_Seq_Analysis";

if (!(-d $RNA_Path))
{
    print "$RNA_Path does not exist. Either it was deleted, moved or the scripts required for it have not been run.\n";
    exit;
}

if(-d "$RNA_Path/peds" and -d "$RNA_Path/maps")
{
            print STDERR "It is going to use the default plink ped and map dirs:\n",
            "$RNA_Path/peds\n","$RNA_Path/maps\n";
}
else
{
            print "There are no peds and maps dirs in the dir $RNA_Path\n";
            exit;
} 

chdir "$RNA_Path";

if(!defined $imputation)
{
    $imputation = "$RNA_Path/phased";
    print STDERR "using $imputation as default\n"; 
}

`mkdir -p $imputation` unless(-d "$imputation");
`rm -f $imputation/*`;

`mkdir "$RNA_Path/logs"` unless(-d "$RNA_Path/logs");
`rm -f $RNA_Path/logs/*`;

my $OneKG_Ref_Path = "$TCGA_Pipeline_Dir/Database/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono";

my $Phased_hap = "$imputation";

my $Impute2out = "$Analysispath/$disease_abbr/phased_imputed_raw_out";

#submit_shapeit(path to ALL.integrated_phase1_SHAPEIT_16-06-14.nomono,path to the RNA_Seq_Analysis directory,user defened directory from command line or default directory)
my @shapeit_cmds = $impute_plink->submit_shapeit("$OneKG_Ref_Path","$RNA_Path","$imputation");

mce_map {
      system("$shapeit_cmds[$_]");
                  } 0..$#shapeit_cmds;

#fetch_Chrom_Sizes(reference genome(e.g. hg19))
$impute_plink->fetch_Chrom_Sizes("hg19");

`cat chr_lens_grep_chr | grep -v Un | grep -v random | grep -v hap | grep -v M | grep -v Y | grep chr > chr_lens`;

`cat chr_lens|head -n 11 > file_for_submit`;
#submit first 11 chrs for imputation
#submit_all(file with chr sizes,path to ALL.integrated_phase1_SHAPEIT_16-06-14.nomono,user defened directory from command line or default directory,path to phased_imputed_raw_out)
my @imput2cmds = $impute_plink->submit_all("file_for_submit", $OneKG_Ref_Path, $Phased_hap, $Impute2out);

mce_map {
      system("$imput2cmds[$_]");
                  } 0..$#imput2cmds;

#To save disk space;
#remove them permentately!
`rm -f $Impute2out/*_allele_probs`;
`rm -f $Impute2out/*_info_by_sample`;
`rm -f $Impute2out/*_info`;
`rm -f $Impute2out/*_summary`;
`rm -f $Impute2out/*_warnings`;

`cat chr_lens|tail -n 12 > file_for_submit`;
#submit next 12 chrs for imputation
undef @imput2cmds;
@imput2cmds = $impute_plink->submit_all( "file_for_submit", $OneKG_Ref_Path, $Phased_hap, $Impute2out);
mce_map
{
            system("$imput2cmds[$_]");

} 0..$#imput2cmds;

#To save disk space;
#remove them permentately!
`rm -f $Impute2out/*_allele_probs`;
`rm -f $Impute2out/*_info_by_sample`;
`rm -f $Impute2out/*_info`;
`rm -f $Impute2out/*_summary`;
`rm -f $Impute2out/*_warnings`;

print "All jobs are done for $disease_abbr.\n";

$time = localtime;
print "Script finished on $time.\n";

exit;