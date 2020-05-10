# MicrobeAnnotator
Easy-to-use pipeline for the comprehensive metabolic annotation of microbial genomes.

## Table of Content
* [Features] (#features)
* [Citation] (#citation)
* [Requirements] (#requirements)
* [Installation] (#installation)
* [Usage] (#usage)
* [FAQs] (#faqs)
* [License]

- [MicrobeAnnotator](#microbeannotator)
  * [Table of Content](#table-of-content)
  * [Features](#features)
  * [Citation](#citation)


## Features
MicrobeAnnotator uses an iterative approach to annotate microbial genomes (Bacteria, Archaea and Virus) starting from proteins predicted using your favorite ORF prediction tool, e.g. Prodigal.
The iterative approach is composed of three or five main steps (depending on the flavor of MicrobeAnnotator you run).
- Search of initial protein dataset using KOFamscan, extraction of unannotated proteins.
- Search of proteins extracted in the previous step using UniProt Swissprot database. Extract unannotated proteins and repeat search using RefSeq and Trembl (if running the full version).
- Summarize the metabolic potential using KEGG modules by extracting KO numbers associated with each match in the databases used.

## Citation
Comming soon.

## Requirements:
- Aspera Connect\
    To install aspera in a Linux system follow (example with version 3.9.8.176272):\
    wget https://download.asperasoft.com/download/sw/connect/3.9.8/ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.tar.gz\
    tar xvfz ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.tar.gz\
    bash ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.sh\
    (Usually it will be installed in your home under your user in "/home/[user]/.aspera/connect/bin", it will show you where)\
    Take note of this installation folder and, for your peace of mind lets add this folder to your PATH:\
    Add the following line to your ~/.bashrc file: export PATH="$PATH:/home/cruizperez/.aspera/connect/bin"\
    Now reload the ~/.bashrc file with: source ~/.bashrc file
    Now you have installed Aspera Connect!\
- HMMER >= 3.1
- Ruby >= 2.4
- GNU Parallel
  
  
