#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Download the latest uniprot protein datasets from Uniprot and
# the accompanying .dat files necessary to create the database.
# ==============================================================================
"""


# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.logging import setup_logger

from shutil import copyfileobj
from pathlib import Path
from sys import argv

import urllib.request
import argparse
import wget
import gzip
# ==============================================================================


# ==============================================================================
# Initalize logger
# ==============================================================================
logger = setup_logger(__name__)
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
# Function to download Swissprot protein sequences
def swissprot_fasta_downloader(database_directory: Path) -> Path:
    """
    Downloads protein sequences from UniProt SwissProt.

    Args:
        database_directory (Path): Database directory to store files. 

    Returns:
        Path: Fasta file with SwissProt protein sequences.
    """
    # Create protein folder
    merged_db_folder = database_directory / "protein_db"
    merged_db_folder.mkdir(parents=True, exist_ok=True)
    # Get release number
    rel_url = (
        f"ftp://ftp.uniprot.org/pub/databases/uniprot/"
        f"current_release/knowledgebase/complete/reldate.txt")
    release_file = urllib.request.urlopen(rel_url)
    release_number = release_file.readline().decode('utf-8').split()[3]
    # Write release information into file
    with open(merged_db_folder / "uniprot_release.txt", 'w') as relinfo:
        relinfo.write("Release number (date) {}".format(release_number))
    # Download fasta files
    logger.info("Downloading SwissProt Fasta files")
    # Overwrite existing temporal files
    output_fasta = merged_db_folder / "uniprot_sprot.fasta.gz"
    if output_fasta.is_file():
        output_fasta.unlink()
    # Download swissprot
    fasta_url = (
        f"ftp://ftp.uniprot.org/pub/databases/uniprot/"
        f"current_release/knowledgebase/complete/uniprot_sprot.fasta.gz")
    wget.download(fasta_url, out=str(output_fasta))
    decompressed_fasta = merged_db_folder / "uniprot_sprot.fasta"
    with gzip.open(output_fasta, 'rt') as compressed_fh, \
        open(decompressed_fasta, 'w') as decompressed_fh:
        copyfileobj(compressed_fh, decompressed_fh)
        output_fasta.unlink()
    logger.info("Finished")

    return output_fasta

# Function to download SwissProt annotations
def swissprot_dat_downloader(database_directory: Path) -> Path:
    """
    Downloads UniProt SwissProt protein annotations in the form of a .dat file.

    Args:
        database_directory (Path): Database directory to store files.

    Returns:
        Path: Path to downloaded .dat annotation file.
    """
    # Create dat files folder
    temp_dat_files = database_directory / "temp_swissprot_dat_files"
    temp_dat_files.mkdir(parents=True, exist_ok=True)
    # Overwrite existing temporal files
    output_dat_file = temp_dat_files / "uniprot_sprot.dat.gz"
    if output_dat_file.is_file():
        output_dat_file.unlink()
    # Download swissprot and trembl dat files
    logger.info("Downloading Swissprot .dat files")
    dat_url = (
        f"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        f"knowledgebase/complete/uniprot_sprot.dat.gz"
    )
    wget.download(dat_url, out=str(output_dat_file))
    logger.info("Finished")

    return output_dat_file

# Function to download TrEMBL protein sequences
def trembl_fasta_downloader(database_directory: Path) -> Path:
    """
    Downloads protein sequences from UniProt TrEMBL.

    Args:
        database_directory (Path): Database directory to store files. 

    Returns:
        Path: Fasta file with TrEMBL protein sequences.
    """
    # Create protein folder
    merged_db_folder = database_directory / "protein_db"
    merged_db_folder.mkdir(parents=True, exist_ok=True)
    # Get release number
    rel_url = (
        f"ftp://ftp.uniprot.org/pub/databases/uniprot/"
        f"current_release/knowledgebase/complete/reldate.txt")
    release_file = urllib.request.urlopen(rel_url)
    release_number = release_file.readline().decode('utf-8').split()[3]
    # Write release information into file
    with open(merged_db_folder / "uniprot_release.txt", 'w') as relinfo:
        relinfo.write("Release number (date) {}".format(release_number))
    # Download fasta files
    logger.info("Downloading TrEMBL Fasta files")
    # Overwrite existing temporal files
    output_fasta = merged_db_folder / "uniprot_trembl.fasta.gz"
    if output_fasta.is_file():
        output_fasta.unlink()
    # Download swissprot
    fasta_url = (
        f"ftp://ftp.uniprot.org/pub/databases/uniprot/"
        f"current_release/knowledgebase/complete/uniprot_trembl.fasta.gz")
    wget.download(fasta_url, out=str(output_fasta))
    decompressed_fasta = merged_db_folder / "uniprot_trembl.fasta"
    with gzip.open(output_fasta, 'rt') as compressed_fh, \
        open(decompressed_fasta, 'w') as decompressed_fh:
        copyfileobj(compressed_fh, decompressed_fh)
        output_fasta.unlink()
    logger.info("Finished")

# Function to download TrEMBL annotations
def trembl_dat_downloader(database_directory: Path) -> Path:
    """
    Downloads UniProt TrEMBL protein annotations in the form of a .dat file.

    Args:
        database_directory (Path): Database directory to store files.

    Returns:
        Path: Path to downloaded .dat annotation file.
    """
    # Create dat files folder
    temp_dat_files = database_directory / "temp_trembl_dat_files"
    temp_dat_files.mkdir(parents=True, exist_ok=True)
    # Overwrite existing temporal files
    output_dat_file = temp_dat_files / "uniprot_trembl.dat.gz"
    if output_dat_file.is_file():
        output_dat_file.unlink()
    # Download swissprot and trembl dat files
    logger.info("Downloading TrEMBL .dat files")
    dat_url = (
        f"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/"
        f"knowledgebase/complete/uniprot_trembl.dat.gz"
    )
    wget.download(dat_url, out=str(output_dat_file))
    logger.info("Finished")

    return output_dat_file
# ==============================================================================


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            f"This script downloads the Fasta and UniProt.dat files\n"
            f"that are needed to build the Uniprot annotation database.\n"
            f"By default it downloads both but either can be specified.\n"
            f"Mandatory parameters: -o [output directory]\n"
            f"Optional parameters: See {argv[0]} -h\n"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory")
    mandatory_arguments.add_argument(
        '-o', '--outdir', dest='outdir',
        action='store', required=True,
        help='Output directory to store the fasta and .dat files.')
    # Setup optional arguments
    optional_arguments = parser.add_argument_group("Optional")
    optional_arguments.add_argument(
        '--light', dest='light', action='store_true', required=False,
        help='Light mode, only downloads the swissprot-related data.')
    optional_arguments.add_argument(
        '--proteins', dest='proteins', action='store_true', required=False,
        help='Only download the protein Fasta files.')
    optional_arguments.add_argument(
        '--datfile', dest='datfile', action='store_true', required=False,
        help='Only download the .dat files.')
    # If no arguments are provided
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    # Parse arguments
    outdir = arguments.outdir
    proteins = arguments.proteins
    datfile = arguments.datfile
    light = arguments.light

    # Run functions
    if light:
        if proteins:
            swissprot_fasta_downloader(outdir)
        elif datfile:
            swissprot_dat_downloader(outdir)
        else:
            swissprot_fasta_downloader(outdir)
            swissprot_dat_downloader(outdir)
    else:
        if proteins:
            swissprot_fasta_downloader(outdir)
            trembl_fasta_downloader(outdir)
        elif datfile:
            swissprot_dat_downloader(outdir)
            trembl_dat_downloader(outdir)
        else:
            swissprot_fasta_downloader(outdir)
            trembl_fasta_downloader(outdir)
            swissprot_dat_downloader(outdir)
            trembl_dat_downloader(outdir)
# ==============================================================================


# ==============================================================================
# Run main function
# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
