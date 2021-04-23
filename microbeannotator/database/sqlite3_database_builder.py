#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Create SQLite database from a tab-delimited tables.
# Used to create annotation parsing databases.
# ==============================================================================
"""

# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.logging import setup_logger

from pathlib import Path
from sys import exit
from sys import argv

import argparse
import sqlite3
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
# Function to create SwissProt table
def create_swissprot_table(database: str, final_table: Path) -> None:
    """
    Creates and adds SwissProt table to SQLite database.

    Args:
        database (str): Path to SQLite database.
        final_table (Path): Table to read and add to database.
    """
    logger.info("Creating RefSeq table")
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS swissprot")
    cursor.execute('CREATE TABLE swissprot \
        (gene_id TEXT, accession TEXT, gene_name TEXT, ko_number TEXT, \
        organism TEXT, taxonomy TEXT, function TEXT, compartment TEXT, \
        process TEXT, interpro TEXT, pfam TEXT, ec_number TEXT)')
    with open(final_table, 'r') as swissprot:
        record_counter = 0
        records = []
        for line in swissprot:
            # Commit changes after 500000 records
            if record_counter == 500000:
                cursor.execute("begin")
                cursor.executemany('INSERT INTO swissprot VALUES(?, ?, ?, ?, \
                    ?, ?, ?, ?, ?, ?, ?, ?)', records)
                cursor.execute("commit")
                record_counter = 0
                records = []
            else:
                records.append(tuple(line.rstrip('\n').split("\t")))
                record_counter += 1
        # Commit remaining records
        if record_counter > 0:
            cursor.execute("begin")
            cursor.executemany('INSERT INTO swissprot VALUES(?, ?, ?, ?, ?, \
                ?, ?, ?, ?, ?, ?, ?)', records)
            cursor.execute("commit")
        # Create index for faster access
        cursor.execute('CREATE INDEX swissprot_id ON swissprot (gene_id)')
    logger.info("Finished")

def create_trembl_table(database: str, final_table: Path) -> None:
    """
    Creates and adds TrEMBL table to SQLite database.

    Args:
        database (str): Path to SQLite database.
        final_table (Path): Table to read and add to database.
    """
    logger.info("Creating TrEMBL table")
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS trembl")
    cursor.execute('CREATE TABLE trembl \
        (gene_id TEXT, accession TEXT, gene_name TEXT, ko_number TEXT, \
        organism TEXT, taxonomy TEXT, function TEXT, compartment TEXT, \
        process TEXT, interpro TEXT, pfam TEXT, ec_number)')
    with open(final_table, 'r') as trembl:
        record_counter = 0
        records = []
        for line in trembl:
            # Commit changes after 500000 records
            if record_counter == 500000:
                cursor.execute("begin")
                cursor.executemany('INSERT INTO trembl VALUES(?, ?, ?, ?, ?, \
                    ?, ?, ?, ?, ?, ?, ?)', records)
                cursor.execute("commit")
                record_counter = 0
                records = []
            else:
                records.append(tuple(line.rstrip('\n').split("\t")))
                record_counter += 1
        # Commit remaining records
        if record_counter > 0:
            cursor.execute("begin")
            cursor.executemany('INSERT INTO trembl VALUES(?, ?, ?, ?, ?, \
                ?, ?, ?, ?, ?, ?, ?)', records)
            cursor.execute("commit")
        # Create index for faster access
        cursor.execute('CREATE INDEX trembl_id ON trembl (gene_id)')
    logger.info("Finished")

def create_refseq_table(database: str, final_table: Path) -> None:
    """
    Creates and adds RefSeq table to SQLite database.

    Args:
        database (str): Path to SQLite database.
        final_table (Path): Table to read and add to database.
    """
    # Store protein ko annotations
    # Download refseq KO annotations
    logger.info("Parsing RefSeq KO numbers")
    annot_url = (
        f"http://enve-omics.ce.gatech.edu/data/public_microbeannotator/"
        f"01.RefSeq_KEGG_Annotations.txt.gz")
    annotation_file = wget.download(
        annot_url, out="refseq_ko_annotations.tar.gz")
    protein_annotations = {}
    with gzip.open(annotation_file, 'rt') as annotations:
        for line in annotations:
            line = line.rstrip('\n').split("\t")
            protein_annotations[line[0]] = [line[1], line[2]]
    Path(annotation_file).unlink()
    logger.info("Finished")
    # Connect to db and create table
    logger.info("Creating RefSeq table")
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS refseq")
    cursor.execute('CREATE TABLE refseq \
        (gene_id TEXT, product TEXT, taxonomy TEXT, \
        ec_number TEXT, ko_number TEXT, ko_product TEXT)')
    with open(final_table, 'r') as refseq:
        record_counter = 0
        records = []
        for line in refseq:
            # Commit changes after 500000 records
            if record_counter == 500000:
                cursor.execute("begin")
                cursor.executemany(
                    'INSERT INTO refseq VALUES(?, ?, ?, ?, ?, ?)', records)
                cursor.execute("commit")
                record_counter = 0
                records = []
            else:
                line = line.rstrip('\n').split("\t")
                if line[0] in protein_annotations:
                    ko_number = protein_annotations[line[0]][0]
                    ko_product = protein_annotations[line[0]][1]
                    records.append((line[0], line[1], line[2], 
                        line[3], ko_number, ko_product))
                    record_counter += 1
                else:
                    records.append((line[0], line[1], line[2],
                        line[3], "NA", "NA"))
                    record_counter += 1
        # Commit remaining records
        if record_counter > 0:
            cursor.execute("begin")
            cursor.executemany(
                'INSERT INTO refseq VALUES(?, ?, ?, ?, ?, ?)', records)
            cursor.execute("commit")
        # Create index for faster access
        cursor.execute('CREATE INDEX refseq_id ON refseq (gene_id)')
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
            f"This script adds tab-separated tables to a SQLite database.\n"
            f"By default it assumes the first line of the input table has\n"
            f"headers. If not present, they must be provided using the\n"
            f"-- header flag as --header header1,header2,header3, etc."
            f"Mandatory parameters: -i [input table] -d [sqlite database]\n"
            f"-t [table type]\n"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory")
    mandatory_arguments.add_argument(
        '-i', '--input', dest='input_table', action='store', required=True,
        help='Input tab-delimited table to parse, by default assumes headers are present')
    mandatory_arguments.add_argument(
        '-d', '--database', dest='database', action='store', required=True,
        help='Database name, can be exisiting or new.')
    mandatory_arguments.add_argument(
        '-t', '--type', dest='db_type', action='store', required=True,
        help='Table type, can be one of swissprot, trembl, or refseq.')
    # If no arguments are provided
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    # Parse arguments
    input_table = arguments.input_table
    database = arguments.database
    db_type = arguments.db_type
    db_type = db_type.lower()
    # Validate requested type
    type_choice = ['swissprot', 'trembl', 'refseq']
    if db_type not in type_choice:
        logger.error(f"Unvalid database type requested, {db_type}.")
        exit(1)
    # Run functions
    logger.info("Building database table")
    if db_type == "swissprot":
        create_swissprot_table(database, input_table)
    elif db_type == "trembl":
        create_trembl_table(database, input_table)
    elif db_type == "refseq":
        create_refseq_table(database, input_table)
    else:
        exit("Please specify a table type, swissprot, trembl or refseq")
# ==============================================================================


# ==============================================================================
# Run main function
# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
