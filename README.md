# Driver-ASE

A pipeline that correlates somatic mutation in tumor/normal WGS paired samples with gene-level allele-specific expression in an effort to systematically uncover cis-regulatory variants involved in tumorigenesis.

# Introduction

It is well known that some DNA mutations can induce uncontrolled cell division and ultimately lead to cancer. Advances in sequencing over the last five years have enabled researchers to scan all of the DNA within a tumor to identify these mutations, and indeed we are finding thousands, leaving us with the daunting task of distinguishing the potentially few causal changes from the vast majority of “passengers”. Researchers are focusing on protein-coding regions where functional effects can be readily predicted. Unfortunately, this approach misses the expansive oncogenic effects caused by changes in gene regulation (turning genes up/down), which are much more difficult to identify. We have pioneered a new bioinformatic approach to probe high-throughput RNA sequencing (RNA-Seq) data, generated from the molecules that carry out gene expression, to identify the truly causal regulatory mutations. We applied our method to all suitable breast cancer samples in The Cancer Genome Atlas (TCGA) and mapped over a hundred genes whose regulation is perturbed by regulatory mutations in breast tumors.

The pipeline performs ASE analysis with using RNA-Seq, affymetrix genotyping array, whole genome sequence, copy number variant, and methylation data from TCGA. It is divided into 5 parts that have 2 or more Perl scripts for each part. They are used to download, run imputation, perform ASE analysis as well as perform WGS analysis. These scripts are easy to use and provide options and usages. To speed up the analysis, the perl module MCE (Many-Core Engine) was implemented into the pipeline. Additionally, this pipeline uses other software to process data, including samtools, impute2, shapeit, plink, overlapSelect and VarScan with links found below for each software. 

# Install

This pipeline was coded on a Centos 7 x86_64 system and the some of the package names listed below are specific to Centos. If you wish to use the pipeline on a different distrobution, install the equivalent packages for the distrobution that you are using. If you don't want to configure your system to use the pipeline, you could use an image that was made with Docker that has all of the dependencies already installed for the pipeline to run.

Docker Images: https://hub.docker.com/r/mikegreen24/driver-ase/

#### Packages

##### This pipeline also requires some packages to be installed as the other software it uses within it needs them. Some of these packages may be installed already.

openssl - sudo yum -y install openssl

openssl-devel - sudo yum -y install openssl-devel

EPEL repository - sudo yum -y install epel-release

yaml - sudo yum -y install yaml-cpp.x86_64 yaml-cpp-devel.x86_64

gcc - sudo yum -y install gcc

gcc-c++ - sudo yum -y install gcc-c++.x86_64

ncurses-devel - sudo yum -y install ncurses-devel.x86_64

java - yum -y install java

bzip2-devel - yum -y install bzip2-devel

xz-devel - yum -y install xz-devel

libcurl-devel - yum -y install libcurl-devel

perl-core - yum -y install perl-core

firewalld - yum -y install firewalld

curl - yum -y install curl

Developement Tools - yum -y groupinstall "Development Tools"

#### Firewall

Make sure that the firewall allows for incoming and outgoing of traffic for HTTPS/HTTP and FTP.

##### Open these ports to allow communication in an out

systemctl start firewalld

systemctl enable firewalld

firewall-cmd --permanent --add-port=20-21/tcp

firewall-cmd --permanent --add-port=80/tcp

firewall-cmd --permanent --add-port=443/tcp

firewall-cmd --permanent --add-port=3306/tcp

firewall-cmd --reload

Make sure you have these perl modules installed before you run these scripts.

##### In order to install these modules, run the commands below.

Some of the module below may already be installed.

curl -L http://cpanmin.us | perl - YAML

curl -L http://cpanmin.us | perl - FindBin

curl -L http://cpanmin.us | perl - File::Copy

curl -L http://cpanmin.us | perl - File::Basename

curl -L http://cpanmin.us | perl - FileHandle

curl -L http://cpanmin.us | perl - MCE

curl -L http://cpanmin.us | perl - Cwd

curl -L http://cpanmin.us | perl - Math::CDF

curl -L http://cpanmin.us | perl - Getopt::Long

curl -L http://cpanmin.us | perl - LWP::Simple

curl -L http://cpanmin.us | perl - autodie

curl -L http://cpanmin.us | perl - lib

curl -L http://cpanmin.us | perl - strict

curl -L http://cpanmin.us | perl - warnings

##### This pipeline uses curl by default but aria2 is supported.

aria2 - https://aria2.github.io/

##### The commands to install aria2 are as follows:

./configure --with-openssl

sudo make

sudo make install

#### You will need to install the software below as well.

#### NOTE: The URLs provided below may change or the page they direct to may no longer be available. If that happens, the software will need to be aquired elsewhere. Also, the commands below are optional other than the software that requires configure and make to be run.

##### If you don't have root permissions to install the software, contact your system admin to have the software installed. 

#### samtools

download link - http://www.htslib.org/download/

#### Installing samtools

tar -jxvf samtools-1.3.1.tar.bz2(your version number may be different)

remove or backup the archive.

cd into the extracted archive and run the commands below.

./configure

make

sudo make install

find where samtools was installed.

which samtools


#### shapeit(GLIBC v2.12 if using Centos or other distrobutions that don't support GLIBC version 2.20)

download link - https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html

#### Installing shapeit

tar -zxvf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz

Locate the bin directory that was extracted and cd into it.

Check if shapeit is executable.

ls -l shapeit

If shapeit is not executable, change the permissions.

chmod a+x shapeit

place the shapeit executable in the same location as samtools (might be /usr/bin or /usr/local/bin/).

remove or backup the archive.

#### impute2(the static verison was used for this pipeline) 

download link - https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download

#### Installing impute2

tar -zxvf impute_v2.3.2_x86_64_static.tgz

cd into the extracted archive.

Check if impute2 is executable.

ls -l impute2

Change it if it is not.

chmod a+x impute2

move impute2 to the same directory as samtools (might be /usr/bin or /usr/local/bin/).

rename the directory impute2.

remove or backup the archive.


#### plink1.9 or later (Developement) 

download link - https://www.cog-genomics.org/plink2

#### Installing plink

unzip plink_linux_x86_64.zip

Check if executable.

ls -l plink prettify

change if not.

chmod a+x plink prettify

remove or backup the archive.

Place plink and prettify in the same directory as samtools (might be /usr/bin or /usr/local/bin/).


#### overlapSelect 

download link - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/

Find overlapSelect and download it. After it downloads, change the permissions on it to be executable

chmod a+x overlapSelect

move it to the same directory as samtools (might be /usr/bin or /usr/local/bin/).

#### VarScan(need java installed in order to run it) 

download link - https://github.com/dkoboldt/varscan

place the downloaded software in the same directory as samtools (might be /usr/bin or /usr/local/bin/).

When running VarScan, the path to it will need to specified.

Try running samtools or any other software, if it shows the usage or nothing, then it can be used. If it gives any kind of error, you will need to add the directory path they were installed in to your path or they were not properly installed and don't have the proper permissions.

# Download Package

##### To get the pipeline this way, git needs to be installed. The package will be donwnloaded to the current directory the user is in when they run the command.

git clone https://github.com/MichealRollins-Green/Driver-ASE.git


# After Pipeline is Downloaded

##### If the package was downloaded through github and not through git in the command line interface, follow the two steps below.

#### 1. Copy it to a directory and extract it.
#### 2. Enter the directory and copy the Driver-ASE directory to a new location.

Move the Driver-ASE directory anywhere on the system.

This pipeline also uses many different files which are all in the same directory called Database.
The directory with the files can be downloaded from the link below.

#### Database directory download link

https://mega.nz/#!XFFCTTJR!RSeaBGz60ITDXX7pugovkiXutSYOYmgNhaOaiVy2d1Q

#### Once downloaded place the directory in the Driver-ASE directory and decompress it.

#### GDC Key Required.

You will also need to have a GDC account to get the key for downloading as it does expire over a certain period of time. Just go to: https://gdc-portal.nci.nih.gov/ 

#### 1. Go to login. 
#### 2. Click on your account name in the top right to open the drop down and click download token. 
#### 3. Rename the downloaded token to gdc.key.
#### 4. Place it in the Database directory of the pipeline. When running the download scripts, the path to the key will need to be entered.

This pipeline can be used right away in any part of the system once all of the above has been done. You can also add the path to the sctipts if you wish so you don't have to be in the directory or specify the path to them in order to execute them.

# Scripts

Below are the scripts that are included in this pipeline:

0_Download_SNPArray_From_GDC.pl – Downloads SNP-Array data from gdc. The data types that are used in this pipeline are Genotypes and Copy number estimate.

0_Download_SNPArray_Table_From_GDC.pl - Downloads the SNP-Array data table of a specified cancer type from GDC. Only run this script if the table file is missing from the tables directory of the cancer type.

1.0_Prep_SNPs_for_Imputation.pl – This prepares the files that will be used in the imputation process.

1.1_Birdseed_to_ped_and_maps.pl – Takes the birdseed files that were downloaded from GDC and makes the map and ped files that will be used for imputation.

1.2_Shapeit_and_Imputation.pl – Runs shapeit processes as well as run imputation on the files made in the previous script.

1.3_Make_het_cds.pl – Keeps het SNPS and makes cds sorted files.

2.0_Prep_For_Bad_SNPs_CNVs.pl – Performs lookups on the cds sorted files and gets the normal samples out of them.

2.1_Get_Bad_SNPs_and_CNVs.pl – Gets the bad SNPs and CNVs of the cancer type.

3.0_Download_RNASeq_WGS_and_do_Mpileup.pl – Downloads either RNA-Seq or WGS bams for a cancer type. By default this script is set to download a certain number of bams unless the user specifies the amount that they want to download in the command line, deletes them afterwards and downloads more when done.

3.0_Dwnld_RNASeq_WGS_Table_From_GDC.pl  - Downloads the RNA-Seq or WGS table of specified cancer types from GDC. Only run this script if the table file is missing from the tables directory of the cancer type.

3.1_Gene_Level_ASE_Analysis.pl – Performs gene level ASE Analysis on the samples that had mpileups done on them.

3.2_Export_ASE_Data.pl – Performs the final touches by compressing what will be needed for further analysis and moving them into a specific directory.

4.0_Somatic_Variants.pl – Processes WGS data to get somatic variants.

4.1_Upstream_Downstream_Analysis.pl – Performs upstream and downstream analysis and gets the data ready to be analyzed.

# Modules

Modules that are incorporated in the scripts are listed below:

Parsing_Routines.pm - All scripts use this module.

Dwnld_WGS_RNA.pm - 0_Download_SNPArray_From_GDC.pl, 0_Download_SNPArray_Table_From_GDC.pl, 3.0_Download_RNASeq_WGS_and_do_Mpileup.pl, 3.0_Dwnld_RNASeq_WGS_Table_From_GDC.pl 

Imputation_Plink.pm - 1.0_Prep_TCGA_for_Imputation_and_Plink.pl, 1.1_Birdseed_to_ped_and_maps.pl, 1.2_Shapeit_and_Imputation.pl, 1.3_Make_het_cds.pl

Bad_SNPS_and_CNVS.pm - 2.0_Lookups_for_normal.pl, 2.1_Get_Bad_SNPs_and_CNVs.pl

Gene_Level_ASE_Analysis.pm - 3.1_Gene_Level_ASE_Analysis.pl, 3.2_Export_ASE_data.pl

WGS_Analysis.pm - 4.0_Somatic_Variants.pl, 4.1_Upstream_Downstream_Analysis.pl
