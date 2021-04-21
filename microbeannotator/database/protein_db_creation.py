#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Builds the search databases depending on the search method
# selected. Can be one of blast, sword or diamond. These must be in the PATH
# or the path must be provided.
# ==============================================================================
"""

# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.logging import setup_logger

from pathlib import Path
from shutil import which
from sys import exit
from sys import argv

import subprocess
import argparse
# ==============================================================================


# ==============================================================================
# Initalize logger
# ==============================================================================
logger = setup_logger(__name__)
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
# Function to create blast databases
def blastp_db_creator(input_directory: Path, bin_path: Path = None) -> None:
    logger.info('Building blast databases')
    database_files = []
    if bin_path != None:
        makeblasdb_call = Path(bin_path) / 'makeblastdb'
    else:
        makeblasdb_call = 'makeblastdb'
    if which(makeblasdb_call) == None:
        logger.error(
            f"BLAST's binary {makeblasdb_call} not found.\n"
            f"Plase make sure BLAST is installed and if possible in PATH.\n"
        )
        exit(1)
    else:
        for protein_file in Path(input_directory).iterdir():
            if protein_file.suffix == '.fasta':
                output_file_name = protein_file.with_suffix("")
                database_files.append(output_file_name)
                subprocess.call(
                    [makeblasdb_call, '-in', protein_file, '-dbtype', 'prot',
                    '-out', output_file_name])
    with open(Path(input_directory)/'blast_db.list', 'w') as db_list:
        for file in database_files:
            db_list.write(f"{file.name}\n")
    logger.info('Finished')

# Function to create diamond databases
def diamond_db_creator(input_directory, threads, bin_path=None):
    logger.info('Building diamond databases')
    database_files = []
    if bin_path != None:
        diamond_call = Path(bin_path) / 'diamond'
    else:
        diamond_call = 'diamond'
    if which(diamond_call) == None:
        logger.error(
            f"Diamond's binary {diamond_call} not found.\n"
            f"Plase make sure Diamond is installed and if possible in PATH.\n"
        )
        exit(1)
    else:
        for protein_file in Path(input_directory).iterdir():
            if protein_file.suffix == '.fasta':
                output_file_name = protein_file.with_suffix("")
                final_db_name = Path(output_file_name).with_suffix('.dmnd')
                database_files.append(final_db_name)
                subprocess.call(
                    [diamond_call, 'makedb', '--in', protein_file, '-d',
                    output_file_name, '--threads', str(threads)])
    with open(Path(input_directory)/'diamond_db.list', 'w') as db_list:
        for file in database_files:
            db_list.write(f"{file.name}\n")
    logger.info('Finished')

# Function to create sword database
def sword_db_creator(input_directory, bin_path=None):
    logger.info('Building sword databases')
    logger.info('No DB needed for sword, just checking if sword is in PATH')
    if bin_path != None:
        sword_call = Path(bin_path) / 'sword'
    else:
        sword_call = 'sword'
    if which(sword_call) == None:
        logger.warning(
            f"Sword's binary {sword_call} not found.\n"
            f"Plase make sure Sword is installed and if possible in PATH.\n"
        )
    with open(Path(input_directory)/'sword_db.list', 'w') as db_list:
        for file in Path(input_directory).iterdir():
            if file.suffix == '.fasta':
                db_list.write(f"{file.name}\n")
    logger.info('Finished')
# ==============================================================================


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            f"This script creates the search databases used by tools\n"
            f"in MicrobeAnnotator.\n"
            f"Mandatory parameters: -d [directory proteins] -m [method/tool]\n"
            f"-d [database directory] -s [sqlite database]\n"
            f"Optional parameters: See {argv[0]} -h"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory")
    mandatory_arguments.add_argument(
        '-d', '--dirprot', dest='dirprot', action='store', required=True,
        help='Directory where all raw fasta files are located.')
    mandatory_arguments.add_argument(
        '-m', '--method', dest='method', action='store', required=True,
        help='Search (and DB creation) method. One of blast, diamond or sword')
    # Setup optional arguments
    optional_arguments = parser.add_argument_group("Optional")
    optional_arguments.add_argument(
        '--bin_path', dest='bin_path', action='store', required=False,
        help='Path to binary folder for selected method/tool.')
    # If no arguments are provided
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    # Parse arguments
    dirprot = arguments.dirprot
    method = arguments.method
    bin_path = arguments.bin_path
    method = method.lower()
    # Run functions
    if method == 'blast':
        blastp_db_creator(dirprot, bin_path)
    elif method == 'diamond':
        diamond_db_creator(dirprot, bin_path)
    elif method == 'sword':
        sword_db_creator(dirprot, bin_path)
    else:
        logger.error(
            f"Search method not recognized. Must be blast, diamond, or sword.")
        exit(1)
# ==============================================================================


# ==============================================================================
# Run main function
# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
