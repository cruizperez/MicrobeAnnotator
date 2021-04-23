#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Parses compressed genbank files and extracts protein
information relevant for annotation purposes.
# ==============================================================================
"""

# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.logging import setup_logger

from pathlib import Path
from typing import List
from Bio import SeqIO
from sys import argv

import multiprocessing
import argparse
import shutil
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
# Function to parse Genbank files and create table
def table_creator(genbank_file: Path) -> Path:
    """
    Parses Genbank files and extract annotation information

    Args:
        genbank_file (Path): Path to compressed genbank file

    Returns:
        Path: Path to temporal table created
    """
    temporal_table = Path(genbank_file).with_suffix('.temp')
    with gzip.open(genbank_file, 'rt') as uncompressed_genbank, \
    open(temporal_table, 'w') as output_file:
        for record in SeqIO.parse(uncompressed_genbank, "genbank"):
            protein_id = record.id
            product = record.description.split("[")[0].strip()
            product = product.lower()
            taxonomy = ""
            species = ""
            ec_number = ""
            if 'organism' not in record.annotations:
                species = 'NA'
            else:
                species = record.annotations['organism']
                if 'taxonomy' not in record.annotations and species == "NA":
                    taxonomy = 'NA'
                elif 'taxonomy' not in record.annotations and species != "NA":
                    taxonomy = species
                else:
                    taxonomy = f"{', '.join(record.annotations['taxonomy'])}, {species}" 
                if product == "" or "unknown" in product:
                    for feature in record.features:
                        if feature.type == "Protein":
                            if 'product' not in feature.qualifiers:
                                if 'name' in feature.qualifiers:
                                    product = feature.qualifiers['name'][0]
                                else:
                                    product = 'NA'
                            else:
                                product = feature.qualifiers['product'][0]
            for feature in record.features:
                if feature.type == "Protein" and \
                    "EC_number" in feature.qualifiers:
                    ec_number = feature.qualifiers["EC_number"][0]
                else:
                    ec_number = "NA"
            output_file.write(
                f"{protein_id}\t{product}\t{taxonomy}\t{ec_number}\n")

    return temporal_table

# Function to parse Genbank files using multiprocessing
def table_generator_worker(
    genbank_list: List[Path], threads: int) -> List[Path]:
    """
    Worker calling table_creator function

    Args:
        genbank_list (List[Path]): List of compressed GenBank files
        threads (int): Threads to use

    Returns:
        List[Path]: List of temporal tables created
    """
    logger.info("Processing GenBank files")
    temp_table_list = []
    try:
        pool = multiprocessing.Pool(threads)
        temp_table_list = pool.map(table_creator, genbank_list)
    finally:
        pool.close()
        pool.join()
    logger.info("Finished")

    return temp_table_list

# Funtion to merge temporary files
def table_merger(temp_table_list: List[Path], final_table: Path) -> None:
    """
    Merges individual annotation tables into a single RefSeq annotation table

    Args:
        temp_table_list (List[Path]): List of paths to temporary tables created
        final_table (Path): Final table file
    """
    logger.info("Merging RefSeq annotation tables")
    with open(final_table, 'w') as output_file:
        for file in temp_table_list:
            with open(file) as temp:
                shutil.copyfileobj(temp, output_file)
    logger.info("Finished")
# ==============================================================================


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            f"This script parses compressed RefSeq Genbank files,\n"
            f"extracts relevant information and stores in a single table.\n"
            f"Mandatory parameters: -i [input file(s)] -o [output table]\n"
            f"Optional parameters: See {argv[0]} -h\n"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory")
    mandatory_arguments.add_argument(
        '-i', '--input_files', dest='input_files', action='store',
        nargs='+', required=True,
        help=(
            f"Space-separated list of compressed genbank files,\n"
            f"e.g., '$(ls *.gz).'"))
    mandatory_arguments.add_argument(
        '-o', '--output_file', dest='output_file',
        action='store', required=True,
        help="File to store final table with record information.")
    # Setup optional arguments
    optional_arguments = parser.add_argument_group("Optional")
    optional_arguments.add_argument(
        '-t', '--threads', dest='threads', action='store',
        default=1, type=int, required=False,
        help='Number of threads to use, by default 1.')
    # If no arguments are provided
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    # Parse arguments
    input_files = arguments.input_files
    output_file = arguments.output_file
    threads = arguments.threads
    # Run functions
    temp_tables = table_generator_worker(input_files, threads)
    table_merger(temp_tables, output_file)
# ==============================================================================


# ==============================================================================
# Run main function
# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
