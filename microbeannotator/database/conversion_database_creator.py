#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Downloads the information necessary to create interconversion of
# accessions and appends it into the MicrobeAnnotator database.
# ==============================================================================
"""

# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.logging import setup_logger

from pathlib import Path
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
# Create table with links from RefSeq -> UniProt
def create_refseq_to_uniprot(
    input_table: Path, sqlitedb: Path) -> None:
    """[summary]

    Args:
        input_table (Path): [description]
        database_to_append (Path): [description]
    """
    logger.info("Creating RefSeq to UniProt table")
    # Connect to database
    conn = sqlite3.connect(sqlitedb)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS refseq_to_uniprot")
    cursor.execute('CREATE TABLE refseq_to_uniprot \
        (refseq_id TEXT, uniprot_id TEXT)')
    # Read input table and create table within SQLite database
    with open(input_table, 'r') as refseq_conversion:
        record_counter = 0
        records = []
        for line in refseq_conversion:
            # Commit changes after 500000 records
            if record_counter == 500000:
                cursor.execute("begin")
                cursor.executemany(
                    'INSERT INTO refseq_to_uniprot VALUES(?, ?)', records)
                cursor.execute("commit")
                record_counter = 0
                records = []
            else:
                line = line.strip('\n').split("\t")
                records.append((line[1], line[0]))
                record_counter += 1
        # Commit remaining records
        if record_counter > 0:
            cursor.execute("begin")
            cursor.executemany(
                'INSERT INTO refseq_to_uniprot VALUES(?, ?)', records)
            cursor.execute("commit")
        # Create index for faster access
        cursor.execute(
            'CREATE INDEX refseq_index ON refseq_to_uniprot (refseq_id)')
    logger.info("Finished")

# Function to download and create KO to EC table 
def create_ko_to_ec(database_directory: Path, sqlitedb: Path) -> str:
    """
    Downloads and parses KO to EC identifier assignments and
    creates table in conversion.db

    Args:
        database_directory (Path): Directory for MicrobeAnnotator database.
        sqlitedb (Path): Path to SQLite database.

    Returns:
        str: Temporal KO_to_EC_identifiers.txt.gz file path.
    """
    logger.info("Creating KO to EC table")
    # Download information
    download_output = str(database_directory / "KO_to_EC_identifiers.txt.gz")
    down_url = (
        f"http://enve-omics.ce.gatech.edu/data/public_microbeannotator/"
        f"02.KO_to_EC_identifiers.txt.gz"
    )
    wget.download(down_url, out=download_output)
    # Connect with database
    conn = sqlite3.connect(sqlitedb)
    cursor = conn.cursor()
    # Create table with correspondence KO -> EC
    cursor.execute("DROP TABLE IF EXISTS ko_to_ec")
    cursor.execute('CREATE TABLE ko_to_ec \
        (ko_identifier TEXT, ec_identifier TEXT)')
    with gzip.open(download_output, 'rt') as ko_conversion:
        record_counter = 0
        records = []
        for line in ko_conversion:
            # Commit changes after 500000 records
            if record_counter == 50000:
                cursor.execute("begin")
                cursor.executemany('INSERT INTO ko_to_ec VALUES(?, ?)', records)
                cursor.execute("commit")
                record_counter = 0
                records = []
            else:
                line = line.strip('\n').split("\t")
                records.append((line[0], line[1]))
                record_counter += 1
        # Commit remaining records
        if record_counter > 0:
            cursor.execute("begin")
            cursor.executemany('INSERT INTO ko_to_ec VALUES(?, ?)', records)
            cursor.execute("commit")
        # Create index for faster access
    cursor.execute('CREATE INDEX ko_id ON ko_to_ec (ko_identifier)')
    cursor.execute('CREATE INDEX ec_id ON ko_to_ec (ec_identifier)')
    logger.info("Finished")

    return download_output

# Function to download and parse interpro tables
def create_interpro_tables(database_directory: Path, sqlitedb: Path) -> str:
    """
    Downloads and parses InterPro tables for intercoversion of ids.

    Args:
        database_directory (Path): Directory for MicrobeAnnotator database.
        sqlitedb (Path): Path to SQLite database.

    Returns:
        str: Path to temporal interpro_metadata.xml.gz file.
    """
    logger.info("Creating InterPro tables")
    # Download information
    download_output = str(database_directory / "interpro_metadata.xml.gz")
    wget.download("ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz",
                            out=download_output)
    interpro_to_ec = {}
    # Parse file information into dictionaries to merge into the SQLite database
    with gzip.open(download_output, 'rt') as infile:
        interpro_id = ""
        ec_identifier = []
        for line in infile:
            if line.startswith("<interpro id="):
                interpro_id = line.strip().split('"')[1]
            elif 'db_xref db="EC"' in line:
                identifier = line.strip().split('"')[3]
                ec_identifier.append(identifier)
            elif line.startswith("</interpro>"):
                if len(ec_identifier) > 0:
                    interpro_to_ec[interpro_id] = list(set(ec_identifier))
                interpro_id = ""
                ec_identifier = []

    # Connect with the database
    conn = sqlite3.connect(sqlitedb)
    cursor = conn.cursor()
    # Create table with correspondence InterPro -> EC
    cursor.execute("DROP TABLE IF EXISTS interpro_to_ec")
    cursor.execute('CREATE TABLE interpro_to_ec \
        (interpro_id TEXT, ec_identifier TEXT)')
    record_counter = 0
    records = []
    for interproscan, ec_id in interpro_to_ec.items():
        for ec_id_record in ec_id:
            # Commit changes after 5000 records
            if record_counter == 5000:
                cursor.execute("begin")
                cursor.executemany(
                    'INSERT INTO interpro_to_ec VALUES(?, ?)', records)
                cursor.execute("commit")
                record_counter = 0
                records = []
            else:
                records.append((interproscan, ec_id_record))
                record_counter += 1
    # Commit remaining records
    if record_counter > 0:
        cursor.execute("begin")
        cursor.executemany('INSERT INTO interpro_to_ec VALUES(?, ?)', records)
        cursor.execute("commit")
    # Create index for faster access
    cursor.execute(
        'CREATE INDEX interpro_index_ce ON interpro_to_ec (interpro_id)')
    logger.info("Finished")

    return download_output
# ==============================================================================


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            f"This script downloads and parses tables required for\n"
            f"intercoversion of identifiers and stores in a SQLite database.\n"
            f"Mandatory parameters: -i [input refseq-uniprot table] \n"
            f"-d [database directory] -s [sqlite database]\n"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory")
    mandatory_arguments.add_argument(
        '-i', '--input_table', dest='intable', action='store', required=True,
        help='RefSeq to Uniprot table to create correspondence table.')
    mandatory_arguments.add_argument(
        '-d', '--database', dest='database', action='store', required=True,
        help='Database directory where data should be stored.')
    mandatory_arguments.add_argument(
        '-s', '--sqlitedb', dest='sqlitedb', action='store', required=True,
        help='SQLite database to create tables in.')
    # If no arguments are provided
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    # Parse arguments
    intable = arguments.intable
    database = arguments.database
    sqlitedb = arguments.sqlitedb
    # Run functions
    try:
        create_refseq_to_uniprot(intable, sqlitedb)
    except:
        logger.error("Could not create RefSeq to Uniprot table.")
        exit(1)
    try:
        create_ko_to_ec(database, sqlitedb)
    except:
       logger.error("Could not create KO to EC table.")
       exit(1)
    try:
        create_interpro_tables(database, sqlitedb)
    except:
       logger.error("Could not create InterPro tables.")
       exit(1)
# ==============================================================================


# ==============================================================================
# Run main function
# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
