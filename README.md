# Microbe_Annotator
Pipeline for metabolic annotation of microbial genomes

Requirements:
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
  
  
