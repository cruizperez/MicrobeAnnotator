# MicrobeAnnotator
Easy-to-use pipeline for the comprehensive metabolic annotation of microbial genomes.

## Content Table
  * [Features](#features)
  * [Citation](#citation)
  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [FAQs](#faqs)
  * [License](#license)

## Features
MicrobeAnnotator uses an iterative approach to annotate microbial genomes (Bacteria, Archaea and Virus) starting from proteins predicted using your favorite ORF prediction tool, e.g. Prodigal.
The iterative approach is composed of three or five main steps (depending on the flavor of MicrobeAnnotator you run).
- Search of initial protein dataset using KOFamscan, extraction of unannotated proteins.
- Search of proteins extracted in the previous step using UniProt Swissprot database. Extract unannotated proteins and repeat search using RefSeq and Trembl (if running the full version).
- Summarize the metabolic potential using KEGG modules by extracting KO numbers associated with each match in the databases used. The summary output is a matrix with module completion and two plots showing module completeness per genome (see below).

## Citation
Comming soon.

## Requirements:
- Programs:
   - [Aspera Connect](https://downloads.asperasoft.com/connect2/)
   - [KofamScan](https://github.com/takaram/kofam_scan)
   - [HMMER](http://hmmer.org/)>=3.1
   - [Ruby](https://www.ruby-lang.org/en/)>=2.4
   - [GNU Parallel](https://www.gnu.org/software/parallel/)\
Either:
   - [Blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)>=2.2
   - [Diamond](https://github.com/bbuchfink/diamond)>=0.9
   - [Sword](https://github.com/rvaser/sword)>=1.0.4
- Python Modules:
   - matplotlib
   - seaborn >= 0.10.1
   - pandas
   - argparse
   - pathlib
   - shutil
   - subprocess
   - gzip
   - biopython
   - sqlite3
   - urllib

## Installation
It appears we need a bunch of pre-requisites to run MicrobeAnnotator! No worries, their installation using Conda is quite easy. If you don't have Conda, you can install it as follows:
1. Download Anaconda from https://www.anaconda.com/products/individual.
2. Run `bash Anaconda-latest-Linux-x86_64.sh` and follow the installation instructions.
3. Once installed you can run `conda -V`. You should get the version of conda that you installed.

Now, let's add the conda channels required to install the pre-requisites:\
`conda config --add channels conda-forge`\
`conda config --add channels bioconda`

Then, create an environment for MicrobeAnnotator:\
`conda create -n microbeannotator blast hmmer ruby parallel diamond sword seaborn biopython wget`\
And activate it:\
`conda activate microbeannotator`

This should take care of most of the requirements except for Aspera Connect and KofamScan, which are a little more involved. Let's install those.
- Aspera Connect\
    To install aspera in a Linux system follow (example with version 3.9.8.176272):\
    `wget https://download.asperasoft.com/download/sw/connect/3.9.8/ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.tar.gz`\
    `tar xvfz ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.tar.gz`\
    `bash ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.sh`\
    (Usually it will be installed in your home under your user in `"/home/[user]/.aspera/connect/bin"`, it will show you where)\
    Take note of this installation folder and, for your peace of mind lets add this folder to your PATH:\
    Add the following line to your `~/.bashrc` file: `export PATH="$PATH:/home/cruizperez/.aspera/connect/bin"`\
    Now reload the ~/.bashrc file with: `source ~/.bashrc`
    Now you have installed Aspera Connect!\
- KofamScan\
    First, let's create a folder to download KofamScan and the databases and files it needs (make sure you have enough space for this (~6Gb). Assume I am creating the folder in my user home `/home/[user]`:\
    `mkdir kofamscan`\
    `cd kofamscan`\
    `wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz`\
    `wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz`\
    `wget ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan-1.2.0.tar.gz`\
    Decompress and untar:\
    `gunzip ko_list.gz`\
    `tar xvfz profiles.tar.gz`\
    `tar xvfz kofamscan-1.2.0.tar.gz`\
    `cd kofamscan-1.2.0`
    When you decompress and enter the kofamscan-1.2.0 folder you will find a `config-template.yml` file, which is required for KofamScan to find the databases. We need to copy it and modify by adding the correct paths to the database we just downloaded.\
    `cp config-template.yml config.yml`\
    Open with your favorite text editor:\
    `vim config.yml`\
    You will see something like this:\
    `# Path to your KO-HMM database`\
    `# A database can be a .hmm file, a .hal file or a directory in which`\
    `# .hmm files are. Omit the extension if it is .hal or .hmm file`\
    `# profile: /path/to/your/profile/db`

    `# Path to the KO list file`\
    `# ko_list: /path/to/your/kolist/file`

    `# Path to an executable file of hmmsearch`\
    `# You do not have to set this if it is in your $PATH`\
    `# hmmsearch: /usr/local/bin/hmmsearch`

    `# Path to an executable file of GNU parallel`\
    `# You do not have to set this if it is in your $PATH`\
    `# parallel: /usr/local/bin/parallel`

    `# Number of hmmsearch processes to be run parallelly`\
    `cpu: 8`
    
    
 
  
  
