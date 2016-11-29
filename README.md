# TCGA-ASE-Analysis-Pipeline
TAAP(TCGA ASE Analysis Pipeline): A pipeline to analyze ASEs of various cancer types and answer how ASE is involved

#Introduction

As TCGA is an international project and contains a lot of NGS data for different cancer types, it was our goal to create a pipeline to analyze this data. With this pipeline, It was decided that allelic specific expression (ASE) can be a method to explore TCGA data.

The goal of this pipeline is to do TCGA ASE analysis with RNA-Seq, affymetrix genotyping array, and whole genome sequence. It is divided into many Perl scripts that are used to download, run imputation, perform ASE analysis as well as perform WGS analysis. these scripts are easy to use and provide options and usages. To speed up the analysis, the perl module MCE (Multiple Core Engine was implemented into the pipeline. Additionally, this pipeline uses other software to process data, including samtools, impute2, shapeit, plink, overlapselect and VarScan with links found below for each software. Also included with this pipeline are perl modules that include a bunch of subroutines that were once scripts.

#Install

This pipeline was coded on a Centos 7 x86_64 system and some of the software listed below are for an x86_64 linux distribution and some of the software may not be compatible with i686 or other versions.

Make sure you have these perl modules installed before you run these scripts.

####In order to install these modules, you will need to use cpan.

Before installing this module it is a good idea to check if they are not already installed:
Enter this command into the command line.

perl -e 'use enter module name here;'

If the module is installed, then nothing will happen but if it is not, it will display an error saying that it could not be found.

FindBin

File::Copy

FileHandle

MCE

Cwd

Math::CDF

Getopt::Long

autodie

####You will need to install the latest version of aria2 since this pipeline uses it to download RNA-Seq, Genotyping array and WGS data.

aria2 - https://aria2.github.io/

####if openssl is not installed the you will have to install that as well because it is needed for https support for aria2.

sudo yum install openssl openssl-devel

####The commands to install aria2 are as follows

./configure --with-openssl

make

sudo make install

####This pipeline also uses curl to get data from GDC, so make sure you have curl installed on your system.

sudo yum install curl

You will need to install the software below as well
####Create a bin directory in your home directory and install the software below to it

cd ~

mkdir bin

####NOTE: The URLs provided below may change or the page they direct to may no longer be available, so if that happens the software will need to be aquired elsewhere.

samtools - https://sourceforge.net/projects/samtools/files/samtools/

shapeit(v2.12 static was used for this pipeline) - https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html

impute2(the static verison was used for this pipeline) - https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download

plink(rename to plink if it is named different) - https://www.cog-genomics.org/plink2

overlapselect - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/

VarScan(need java installed in order to run it) - https://sourceforge.net/projects/varscan/files/

You will need to add the software to your path in order to run them directly from the command line.

cd ~ 

vi .bash_profile

####press i to insert into the file

PATH=$PATH:$HOME/bin:$HOME/bin/samtools/:$HOME/bin/shapeit/bin/:$HOME/bin/impute2/:/usr/bin/:$HOME/bin/plink

####press the Esc key

:wq

source ~/.bash_profile

This pipeline also uses many different which are all in the same directory.
The directory with the files can be downloaded from the link below.
Once downloaded place these files in the TCGA_ASE_Analysis/TCGA_Pipeline/ directory and decompress it. 

####Database directory download link

https://mega.nz/#!bZN3hSob!B5ybrcH_4frJfFv24sMns2XEzHPO5aQsqyrUq1MpnKc

You will also need to have a GDC account to get the key for downloading as it does expire over a certain period of time. Just go to https://gdc-portal.nci.nih.gov/ an go to login. after logging into your account click on your logging name to open the drop down and click download token. after the token is done downloading, rename it gdc.key and place it in the where ever you wish though it makes sense to put it in the Database directory of the pipeline. When running the download scripts, the path to the key will need to be entered.

This pipeline can be used right away in any part of the system once all of the above has been done. You can also add the script directory to your path if you wish so you don't have to be in the directory or specify the path to them in order to execute them.

#Scripts

Below are the scripts that are included in this pipeline:

0_Download_SNPArray_From_GDC.pl – Downloads SNP-Array data from gdc. The data types that are used in this pipeline are Genotypes and Copy number estimate.

0_Download_SNPArray_Table_From_GDC.pl - Downloads the SNP-Array data table of a specified cancer type from GDC. Only run this script if the table file is missing from the tables directory of the cancer type.

1.0_Prep_SNPs_for_Imputation.pl – This prepares the files that will be used in the imputation process.

1.1_Birdseed_to_ped_and_maps.pl – Takes the birdseed files that were downloaded from GDC and makes the map and ped files that will be used for imputation.

1.2_Shapeit_and_Imputation.pl – Runs shapeit processes as well as run imputation on the files made in the previous script.

1.3_Make_het_cds.pl – Keeps het SNPS and makes cds sorted files.

2.0_Lookups_for_normal.pl – Performs lookups in the cds files and gets the normal samples out of them.

2.1_Get_Bad_SNPs_and_CNVs.pl – Gets the bad SNPs and CNVs of the cancer type.

3.0_Download_RNASeq_WGS_and_do_Mpileup.pl – Downloads either RNA-Seq or WGS bams for a cancer type. By default this script is set to download a certain number of bams unless the user specifies the amount that they want to download in the command line, deletes them afterwards and downloads more when done.

3.0_Dwnld_Table_4_RNASeq_WGS_From_GDC.pl - Downloads the RNA-Seq or WGS table of specified cancer types from GDC. Only run this script if the table file is missing from the tables directory of the cancer type.

3.1_ASE_Analysis.pl – Performs gene level ASE Analysis on the samples that had mpileups done on them.

3.2_Export_ASE_Data.pl – Performs the final touches by compressing what will be needed for further analysis and moving them into a specific directory.

4.0_Somatic_Variants.pl – Processes WGS data to get somatic variants.

4.1_Upstream_Downstream_Analysis.pl – Performs upstream and downstream analysis and gets the data ready to be analyzed.

#Modules

Modules that are incorporated in the scripts are listed below:

Parsing_Routines.pm - All scripts use this module.

Dwnld_WGS_RNA.pm - 0_Download_SNPArray_From_GDC.pl, 0_Download_SNPArray_Table_From_GDC.pl, 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl, 3.0_Dwnld_Table_4_RNASeq_WGS_From_GDC.pl

Imputation_Plink.pm - 1.0_Prep_TCGA_for_Imputation_and_Plink.pl, 1.1_Birdseed_to_ped_and_maps.pl, 1.2_Shapeit_and_Imputation.pl, 1.3_Make_het_cds.pl

Bad_SNPS_and_CNVS.pm - 2.0_Lookups_for_normal.pl, 2.1_Get_Bad_SNPs_and_CNVs.pl

TCGA_ASE_Analysis.pm - 3.1_ASE_Analysis_for_other_Bams.pl, 3.2_Export_ASE_data.pl

WGS_Analysis.pm - 4.0_Somatic_Variants.pl, 4.1_Upstream_Downstream_Analysis.pl
