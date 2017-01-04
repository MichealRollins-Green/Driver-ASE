# TAAP (TCGA ASE Analysis Pipeline)

A pipeline that correlates somatic mutation in tumor/normal WGS paired samples with gene-level allele-specific expression in an effort to systematically uncover cis-regulatory variants involved in tumorigenesis.

#Introduction

It is well known that some DNA mutations can induce uncontrolled cell division and ultimately lead to cancer. Advances in sequencing over the last five years have enabled researchers to scan all of the DNA within a tumor to identify these mutations, and indeed we are finding thousands, leaving us with the daunting task of distinguishing the potentially few causal changes from the vast majority of “passengers”. Researchers are focusing on protein-coding regions where functional effects can be readily predicted. Unfortunately, this approach misses the expansive oncogenic effects caused by changes in gene regulation (turning genes up/down), which are much more difficult to identify. We have pioneered a new bioinformatic approach to probe high-throughput RNA sequencing (RNA-Seq) data, generated from the molecules that carry out gene expression, to identify the truly causal regulatory mutations. We applied our method to all suitable breast cancer samples in The Cancer Genome Atlas (TCGA) and mapped over a hundred genes whose regulation is perturbed by regulatory mutations in breast tumors.

The pipeline performs ASE analysis with using RNA-Seq, affymetrix genotyping array, and whole genome sequence, copy number variant, and methylation data from TCGA. It is divided into 4 main Perl scripts that are used to download, run imputation, perform ASE analysis as well as perform WGS analysis. These scripts are easy to use and provide options and usages. To speed up the analysis, the perl module MCE (Multiple Core Engine) was implemented into the pipeline. Additionally, this pipeline uses other software to process data, including samtools, impute2, shapeit, plink, overlapSelect and VarScan with links found below for each software. 

#Install

This pipeline was coded on a Centos 7 x86_64 system and some of the software listed below are for an x86_64 linux distribution and some of the software may not be compatible with i686 or other versions.

Make sure you have these perl modules installed before you run these scripts.

####In order to install these modules, you will need to use cpan.

If cpan is not installed, you will need to install it as well.

Before installing these modules, it is a good idea to check if they are not already installed:

If the module is installed, then nothing will happen but if it is not, it will display an error saying that it could not be found.

FindBin - perl -e 'use FindBin;'

File::Copy - perl -e 'use File::Copy;'

FileHandle - perl -e 'use FileHandle;'

MCE - perl -e 'use MCE;'

Cwd - perl -e 'use Cwd;'

Math::CDF - perl -e 'use Math::CDF;'

Getopt::Long - perl -e 'use Getopt::Long;'

autodie - perl -e 'use autodie;'

lib -  perl -e 'use lib;'

strict -  perl -e 'use strict'

warnings -  perl -e 'use warnings;'

####You will need to install the latest version of aria2 since this pipeline uses it to download RNA-Seq, Genotyping array and WGS data.

aria2 - https://aria2.github.io/

####if openssl is not installed then you will have to install that as well because aria2 needs it to enable https support.

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

####If you don't have root permissions to install the software, contact you system admin to have the software installed.

####samtools

download link - https://sourceforge.net/projects/samtools/files/latest/download?source=files

Place the compressed archive in the bin directory and decompress it

####Installing samtools

tar -jxvf samtools-1.3.1.tar.bz2

Change the name of the directory to samtools

remove or backup the archive

cd samtools

./configure

make

sudo make install


####shapeit(v2.12 static was used for this pipeline) 

download link - https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html

Place the compressed archive in the bin directory and decompress it

####Installing shapeit

tar -zxvf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz

remove or backup the archive

mkdir shapeit

mv bin example LICENCE shapeit


####impute2(the static verison was used for this pipeline) 

download link - https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download

Place the compressed archive in the bin directory and decompress it

####Installing impute2

tar -zxvf impute_v2.3.2_x86_64_static.tgz 

rename the directory impute2

remove or backup the archive


####plink(Stable) 

download link - https://www.cog-genomics.org/plink2

Place the compressed archive in the bin directory and decompress it

####Installing plink

unzip plink_linux_x86_64.zip

mkdir plink1(temporary name)

mv plink prettify LICENSE toy.map toy.ped plink1

mv plink1 plink


####overlapSelect 

download link - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/overlapSelect

place the downloaded software in the bin directory


####VarScan(need java installed in order to run it) 

download link - https://sourceforge.net/projects/varscan/files/

place the downloaded software in the bin directory


You will need to add the software to your path in order for the pipeline to be able to run them. VarScan is the ony exception to this as you will need to put the full path to it in order for it to run.

cd ~ 

vi .bash_profile

####press i to insert into the file

PATH=$PATH:$HOME/bin:$HOME/bin/samtools/:$HOME/bin/shapeit/bin/:$HOME/bin/impute2/:/usr/bin/:$HOME/bin/plink

####press the Esc key

:wq

source ~/.bash_profile

This pipeline also uses many different files which are all in the same directory.
The directory with the files can be downloaded from the link below.

####Database directory download link

https://mega.nz/#!bZN3hSob!B5ybrcH_4frJfFv24sMns2XEzHPO5aQsqyrUq1MpnKc

####Once downloaded place the directory in the TCGA_ASE_Analysis/TCGA_Pipeline/ directory and decompress it.

####GDC Key Required

You will also need to have a GDC account to get the key for downloading as it does expire over a certain period of time. Just go to https://gdc-portal.nci.nih.gov/ 

####1. Go to login. 
####2. Click on your logging name to open the drop down and click download token. 
####3. Rename the downloaded token to gdc.key
####4. Place it in the Database directory of the pipeline. When running the download scripts, the path to the key will need to be entered.

This pipeline can be used right away in any part of the system once all of the above has been done. You can also add the script directory to your path if you wish so you don't have to be in the directory or specify the path to them in order to execute them.

#Scripts

Below are the scripts that are included in this pipeline:

0_Download_SNPArray_From_GDC.pl – Downloads SNP-Array data from gdc. The data types that are used in this pipeline are Genotypes and Copy number estimate.

0_Download_SNPArray_Table_From_GDC.pl - Downloads the SNP-Array data table of a specified cancer type from GDC. Only run this script if the table file is missing from the tables directory of the cancer type.

1.0_Prep_SNPs_for_Imputation.pl – This prepares the files that will be used in the imputation process.

1.1_Birdseed_to_ped_and_maps.pl – Takes the birdseed files that were downloaded from GDC and makes the map and ped files that will be used for imputation.

1.2_Shapeit_and_Imputation.pl – Runs shapeit processes as well as run imputation on the files made in the previous script.

1.3_Make_het_cds.pl – Keeps het SNPS and makes cds sorted files.

2.0_Lookups_for_normal.pl – Performs lookups on the cds sorted files and gets the normal samples out of them.

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
