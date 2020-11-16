#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Download the latest uniprot protein datasets from Uniprot and
# the accompanying .dat files necessary to create the database.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import wget
from shutil import copyfileobj
from pathlib import Path
import urllib.request
import gzip

################################################################################
"""---1.0 Define Functions---"""
def uniprot_fasta_downloader(output_file_folder, light):
    # Create protein folder
    merged_db_folder = Path(output_file_folder) / "01.Protein_DB"
    Path(merged_db_folder).mkdir(parents=True, exist_ok=True)
    # Get release number
    release_file = urllib.request.urlopen("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt")
    release_number = release_file.readline().decode('utf-8').split()[3]
    # Write release information into file
    with open(merged_db_folder / "uniprot_release.txt", 'w') as relinfo:
        relinfo.write("Release number (date) {}".format(release_number))
    # Download fasta files
    print("Downloading fasta files...", flush=True)
    # Overwrite existing temporal files
    if light == True:
        if (merged_db_folder / "uniprot_sprot.fasta.gz").is_file():
            Path.unlink(merged_db_folder / "uniprot_sprot.fasta.gz")
    else:
        if (merged_db_folder / "uniprot_sprot.fasta.gz").is_file():
            Path.unlink(merged_db_folder / "uniprot_sprot.fasta.gz")
        if (merged_db_folder / "uniprot_trembl.fasta.gz").is_file():
            Path.unlink(merged_db_folder / "uniprot_trembl.fasta.gz")
    # Download swissprot and trembl proteins
    if light == True:
        wget.download("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
                  out=str(merged_db_folder / "uniprot_sprot.fasta.gz"))
        compressed = merged_db_folder / "uniprot_sprot.fasta.gz"
        decompressed = merged_db_folder / "uniprot_sprot.fasta"
        with gzip.open(compressed, 'rt') as compressed_fh, open(decompressed, 'w') as decompressed_fh:
            copyfileobj(compressed_fh, decompressed_fh)
            Path.unlink(Path(merged_db_folder / "uniprot_sprot.fasta.gz"))
    else:
        wget.download("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
                  out=str(merged_db_folder / "uniprot_sprot.fasta.gz"))
        wget.download("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
                  out=str(merged_db_folder / "uniprot_trembl.fasta.gz"))
        compressed = merged_db_folder / "uniprot_sprot.fasta.gz"
        decompressed = merged_db_folder / "uniprot_sprot.fasta"
        with gzip.open(compressed, 'rt') as compressed_fh, open(decompressed, 'w') as decompressed_fh:
            copyfileobj(compressed_fh, decompressed_fh)
            Path.unlink(Path(merged_db_folder / "uniprot_sprot.fasta.gz"))
        compressed = merged_db_folder / "uniprot_trembl.fasta.gz"
        decompressed = merged_db_folder / "uniprot_trembl.fasta"
        with gzip.open(compressed, 'rt') as compressed_fh, open(decompressed, 'w') as decompressed_fh:
            copyfileobj(compressed_fh, decompressed_fh)
            Path.unlink(Path(merged_db_folder / "uniprot_trembl.fasta.gz"))
    print("\nDone")
    
def uniprot_dat_downloader(output_file_folder, light):
    # Create dat files folder
    temp_dat_files = Path(output_file_folder) / "02.temp_dat_files"
    Path(temp_dat_files).mkdir(parents=True, exist_ok=True)
    print("Downloading .dat files...", flush=True)
    # Overwrite existing temporal files
    if light == True:
        if (temp_dat_files / "uniprot_sprot.dat.gz").is_file():
            Path.unlink(temp_dat_files / "uniprot_sprot.dat.gz")
    else:
        if (temp_dat_files / "uniprot_sprot.dat.gz").is_file():
            Path.unlink(temp_dat_files / "uniprot_sprot.dat.gz")
        if (temp_dat_files / "uniprot_trembl.dat.gz").is_file():
            Path.unlink(temp_dat_files / "uniprot_trembl.dat.gz")
    # Download swissprot and trembl dat files
    if light == True:
        wget.download("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
                    out=str(temp_dat_files / "uniprot_sprot.dat.gz"))
    else:
        wget.download("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
                    out=str(temp_dat_files / "uniprot_sprot.dat.gz"))
        wget.download("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz",
                      out=str(temp_dat_files / "uniprot_trembl.dat.gz"))
    print("\nDone")


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script downloads the fasta and .dat files that\n'''
            '''are needed to build the Uniprot annotation database. By default it\n'''
            '''downloads both but you can specify either\n'''
            '''Usage: ''' + sys.argv[0] + ''' -f [output_file folder]\n'''
            '''Global mandatory parameters: -f [output_file folder]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-f', '--folder', dest='folder', action='store', required=False,
                        help='Folder to store the fasta and genbank files')
    parser.add_argument('--light', dest='light', action='store_true', required=False,
                        help='Light and fast version, only downloads the swissprot database')
    parser.add_argument('--proteins', dest='proteins', action='store_true', required=False,
                        help='Only download the protein fasta files.')
    parser.add_argument('--datfile', dest='datfile', action='store_true', required=False,
                        help='Only download the genbank files')
    args = parser.parse_args()

    folder = args.folder
    proteins = args.proteins
    datfile = args.datfile
    light = args.light

    # ----------------------------
    if proteins == False and datfile == False:
        uniprot_fasta_downloader(folder, light)
        uniprot_dat_downloader(folder, light)
    elif proteins == True:
        uniprot_fasta_downloader(folder, light)
    elif datfile == True:
        uniprot_dat_downloader(folder, light)
    # ----------------------------

if __name__ == "__main__":
    main()