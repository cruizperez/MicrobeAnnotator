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
- In the end, the program will check those entries that had an annotation match but no KO number and will look for other database identifiers (E.C. and InterPro) in the annotation metadata and will try to link those to KO numbers to refine the annotations. 
- Summarize the metabolic potential using KEGG modules by extracting KO numbers associated with each match in the databases used. The summary output is a matrix with module completion and two plots showing module completeness per genome (see below).

## Citation
MicrobeAnnotator: a user-friendly, comprehensive microbial genome annotation pipeline.
https://doi.org/10.1101/2020.07.20.211847

## Requirements:
- Programs:
   - [Aspera Connect](https://downloads.asperasoft.com/connect2/)
   - [KofamScan](https://github.com/takaram/kofam_scan)
   - [HMMER](http://hmmer.org/) >= 3.1
   - [Ruby](https://www.ruby-lang.org/en/) >= 2.5
   - [GNU Parallel](https://www.gnu.org/software/parallel/)\
Either (or if you prefer all):
   - [Blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) >= 2.2
   - [Diamond](https://github.com/bbuchfink/diamond) >= 0.9
   - [Sword](https://github.com/rvaser/sword) >= 1.0.4
- Python >=3.6,<3.9
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
   - pywget

## Installation
### Conda Installation
It appears we need a bunch of pre-requisites to run MicrobeAnnotator! No worries, their installation using Conda is quite easy. If you don't have Conda, you can install it as follows:
1. Download Anaconda from https://www.anaconda.com/products/individual.
2. Run `bash Anaconda-latest-Linux-x86_64.sh` and follow the installation instructions.
3. Once installed you can run `conda -V`. You should get the version of conda that you installed.

Now, let's add the conda channels required to install the pre-requisites:

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels cruizperez
```

Then, create an environment for MicrobeAnnotator:

```bash
conda create -n microbeannotator microbeannotator blast hmmer ruby=2.5.1 parallel diamond sword seaborn biopython pywget
```

And activate it:

```bash
conda activate microbeannotator
```

Both main scripts (microbeannotator and microbeannotator_db_builder) should be in your path ready for use!
This should take care of most of the requirements except for Aspera Connect and KofamScan, which are a little more involved. Let's install those.

### Pip Installation
Once you have installed the pre-requisites to run MicrobeAnnotator, or if you already had them and you are not using Conda, you can install MicrobeAnnotator using pip:

```bash
pip install microbeannotator
```
Both main scripts (microbeannotator and microbeannotator_db_builder) should be in your path ready for use!

### Aspera Connect
To install aspera in a Linux system follow (example with version 3.9.8.176272):

```bash
wget https://download.asperasoft.com/download/sw/connect/3.9.8/ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.tar.gz
tar xvfz ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.tar.gz
bash ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.sh
```

(Usually it will be installed in your home under your user in `"/home/[user]/.aspera/connect/bin"`, it will show you where)

Take note of this installation folder and, for your peace of mind lets add this folder to your PATH:

Add the following line to your `~/.bashrc` file: `export PATH="$PATH:/home/[user]/.aspera/connect/bin"`\
Note: Replace [user] by your username

Now reload the ~/.bashrc file with: `source ~/.bashrc`

Now you have installed Aspera Connect!

If you cannot install Aspera Connect for some reason, MicrobeAnnotator will use another method to download your data, no worries!

### KofamScan
First, let's create a folder to download KofamScan and the databases and files it needs (make sure you have enough space for this (~6Gb). Assume I am creating the folder in my user home `/home/[your_user]` (note your_user is a placeholder for your username). In addition, note that the version of kofamscan can change in the original location (for example kofam_scan-1.3.0.tar.gz becomes kofam_scan-1.4.0.tar.gz) and therefore you need to adjust your commands to reflect this:

```bash
mkdir kofamscan
cd kofamscan
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
wget ftp://ftp.genome.jp/pub/tools/kofam_scan/kofam_scan-1.3.0.tar.gz

# Decompress and untar:
gunzip ko_list.gz
tar xvfz profiles.tar.gz
tar xvfz kofamscan-1.3.0.tar.gz
cd kofamscan-1.3.0
```

When you decompress and enter the kofamscan-1.3.0 folder you will find a `config-template.yml` file, which is required for KofamScan to find the databases. We need to copy it and modify by adding the correct paths to the database we just downloaded.

```bash
cp config-template.yml config.yml
```

Open with your favorite text editor: `vim config.yml`

You will see something like this:
```yaml
# Path to your KO-HMM database
# A database can be a .hmm file, a .hal file or a directory in which
# .hmm files are. Omit the extension if it is .hal or .hmm file
# profile: /path/to/your/profile/db

# Path to the KO list file
# ko_list: /path/to/your/kolist/file

# Path to an executable file of hmmsearch
# You do not have to set this if it is in your $PATH
# hmmsearch: /usr/local/bin/hmmsearch

# Path to an executable file of GNU parallel
# You do not have to set this if it is in your $PATH
# parallel: /usr/local/bin/parallel

# Number of hmmsearch processes to be run parallelly
cpu: 8
```

This is the information we need to edit with the path to the folder where we downloaded the data:
In this example we will replace: `# profile: /path/to/your/profile/db` by `profile: /home/[your_user]/kofamscan/profiles/prokaryote`

NOTE: If you need to perform Eukaryote searches, you can change the "prokaryote" part by "eukaryote".

Next, we replace: `# ko_list: /path/to/your/kolist/file` by `ko_list: /home/[your_user]/kofamscan/ko_list`

And, given that we have already installed hmmer and parallel earlier we can ignore the other lines.
Finally, you can either add the KofamScan location folder to your $PATH as before with Aspera Connect or later pass it to MicrobeAnnotator when running.

Now, KofamScan is ready to be used! You can test it by running `./exec_annotation -h`.

### Downloading MicrobeAnnotator
Next, we need to download MicrobeAnnotator:

```bash
git clone https://github.com/cruizperez/MicrobeAnnotator.git
```

This should download the program in the folder you are in and everything should be good to go! :)

## Usage
### Database creation
First things first. There are two "flavors" of MicrobeAnnotator depending on your storage and computational capabilities (and your time); the regular and light versions.
The difference between the two is that the regular implementation will make use of four different databases to annotate your proteins, i.e., the KofamScan database and the Swissprot, Trembl and RefSeq databases. On the other hand, the light version will only use KofamScan and Swissprot, which in most cases will annotate the majority of your proteins and will take way less time and space.\
You can decide of course which version to run at any time, but for the time being let's have an example with the full version (to run the light version just add the `--light` flag when calling the script.

The first step is to download and format the databases we want to use for MicrobeAnnotator. For this we will execute the `microbeannotator_db_builder` script within the MicrobeAnnotator folder. You can add see all options and inputs for the script with  `microbeannotator_db_builder -h`.\
Run the script as:\
`microbeannotator_db_builder -d MicrobeAnnotator_DB -m [blast,diamond,sword] -t [# threads]`\
The options we gave the script were:\
`-d MicrobeAnnotator_DB` that is the folder in which all databases will be downloaded and stored.\
`-m [blast,diamond,sword]` the search method you intend to use. For instance, if you choose "blast" the program will format you database to be searched using NCBI blastp. If you noticed we already installed all three options so you can chooce either (note that sword has some specific processor requirements and might not work in older computers, but you can test it before choosing).\
`-t [# threads]` will be the number of processors used to format the databases.\
`--step [#]` will tell MicrobeAnnotator which step to start from (for example if your run failed for some reason). The steps available are: 1. Download data, 2. Parse annotation data, 3. Building SQLite DB, 4. Create interconversion tables, and 5. Build protein DBs.\
`--no_aspera` will tell MicrobeAnnotator that you don't have Aspera Connect installed so it will use another download method.\
`--keep` will tell MicrobeAnnotator to not remove intermediary files (not recommended because it's not necessary and it will take more of your disk space.

If you run the full version (no `--light` flag), you will need at least approximately ~230Gb of space. This can increase depending on the search method used (blast and diamond require the raw fasta databases to be further formated).
If everything went right, you should find inside the MicrobeAnnotator_DB folder, a directory `01.Protein_DB` and a file `02.MicrobeAnnotator.db` that contain all information required by the program to search your proteins against these databases.
If anything went wrong, or if you changed your mind about your searching method after it finished, you can restart the process at any step with the `--step` flag, the help will tell you which are the steps.

### Annotation
This is the main part of MicrobeAnnotator. Here you will search the proteins you provide against the databases we created before and you will receive an annotation table per protein file along with the individual results per database used (in case you want to further explore those), a matrix summarizing the KEGG module completeness and two plots showing said module completeness for an easy, visual comparison between genomes (or files).\
This implies that you can pass multiple protein fasta files to MicrobeAnnotator and they will all be processed in a single script! Awesome right?\
For this annotation step we will use the `MicrobeAnnotator` script (see `/path/to/MicrobeAnnotator -h` for all options).\
Run the script as:\
`microbeannotator -i [fasta_1.fa fasta_2.fa] -d /path/to/MicrobeAnnotator -o [output folder] -m [blast,diamond,sword] -p [# processes] -t [# threads]`\
The options we gave the script were:\
`-i [fasta_1.fa fasta_2.fa]`, the input fasta files (separated by spaces). If you have many, you can pass them as `-i $(ls *.fasta)`.\
`-l [file]`, a file containing the files to process. If you have many files it is easier compared to `-i`\
`-d /path/to/MicrobeAnnotator_DB`, the folder where you created the databases in the previous step.\
`-o [output folder]`, the folder to store all results generated by MicrobeAnnotator.\
`-m [blast,diamond,sword]` the search method you intend to use.\
`-p [# processes]`, refers to the number of protein files to be processed simultaneously, e.g `-p 3` will process three protein files at the same time.\
`-t [# threads]`, refers to the number of processors to use per protein file. For example `-t 5` will use five processors per each protein file.\
Therefore, of you use `-p 2 -t 5`, you will need 10 cores (five for each protein file).\
`--refine`, will tell MicrobeAnnotator to completement the initial annotations using E.C. numbers or InterPro ids and convert them to KO numbers.\
`--continue_run`, we know that MicrobeAnnotator can fail. If you have many files it can be really frustrating having to start over (we have been there :)). That is why we added an option to resume the process from the last completed point. This is the flag that does it. 

At the end you should see your output folder with: \
`annotation_results/` folder: Where you will find the annotation tables per file and the list of KO numbers present in the file.\
`kofamscan_results/` folder: Raw and filtered KofamScan search results.\
`[database]_results/` folder: Raw and filtered search results for each database used (depends of the flavor used).\
A `metabolic_summary_barplot` file: Summarizes the KEGG module completeness per genome (protein file).\
A `metabolic_summary_heatmap` file: Summarizes the KEGG module completeness of all genomes in an-easier-to-understand heatmap (compared to a raw table, that is).\
A `metabolic_summary_module_completeness.tab` file: Matrix describing the completness of all KEGG modules per genome (protein file).

This should be the end result of MicrobeAnnotator, where you can easily compare the metabolic potential of your microbial genomes and if you like, have more information of the annotation of each protein in each file. Happy annotation!

## FAQs
- How long will MicrobeAnnotator take to build the databases?\
    Well, it depends on the version of database you are building. For instance, the full version needs to download the Swissprot, Trembl, and RefSeq databases (and their associated metadata) and parse them into the SQLite database MicrobeAnnotator uses. Depending on your internet connection speed this process can take ~24h (using a single thread), including downloading and processing the data. Fortunately, you just need to do this once! Now, if you are using sword as your search method, you are done, but if you are using Blast or Diamond, there is additional time required for the protein fasta files to be converted into the appropriate formats. This formatting process will take approximately ~5 hours for Blast and ~2.5 hours for Diamond. If you are running the "light" version of the program with sword as your search method, this entire process will only take ~5 minutes (depending of your internet connection speed).
    
- How much space will the database take in my disk?\
    This once again depends on the verison of the program you indend to run and the search method you want to use. For reference the space required by the full database will approximatelly be: 
    - Protein fasta files: ~144Gb
    - Blast-formatted files: ~177Gb
    - Diamond-formatted files: ~148Gb
    - SQLite3 database: ~93Gb\
    In summary, if you build the full database you need at least ~237Gb of free space in your disk. The light version only will take ~0.65Gb (with sword as search method).

- Can I use the light (--light) version of the program to annotate if I built the database using the full version?\
    Of course! You can search and annotate your genomes (protein files) using only kofamscan and the swissprot database, i.e. the light version even if you have the full database available. However, the contrary cannot be done. If you try to perform the full search with the light version of the database, MicrobeAnnotator will not be able to find the required databases and will exit with an error. If you have anough time (and space) you can build the full database and then use the version that better suits your needs at each particular time.

- I have identifiers from other annotation tools, can I use them in MicrobeAnnotator to summarize the metabolic potential?\
    Of course! We realize we are not the only annotation tool out there, so we have an additional script `identifier_conversion.py`, that will allow you to convert ids between databases, either to summarize within microbeannotator or to compare several annotation methods. With this script you can convert between: E.C. <-> KO, InterPro <-> KO, and InterPro <-> E.C. If you have a list of KO identifiers in the end, you can summarize them using another script `ko_mapper.py`.

- I have a list of KOs obtained somewhere. Can I summarize them and skip all the annotation steps?
    Sure! We have made an additional script that does exactly this. You will find it as `ko_mapper.py`.

- I don't see my question here. Where can I ask a question?
    You can submit a ticket in the issues tab and I will do my best to solve it and answer your question. You can also ask for features you would like to see in the future!


## License
See LICENSE


    
    
 
  
  
