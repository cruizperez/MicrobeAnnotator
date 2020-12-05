#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Downloads the information necessary to create interconversion of
# accessions and appends it into the MicrobeAnnotator database.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import wget
from pathlib import Path
import gzip
import sqlite3

################################################################################
"""---1.0 Define Functions---"""
def create_refseq_to_uniprot(input_table, database_to_append, keep):
    # Connect to database
    conn = sqlite3.connect(database_to_append)
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
                cursor.executemany('INSERT INTO refseq_to_uniprot VALUES(?, ?)', records)
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
            cursor.executemany('INSERT INTO refseq_to_uniprot VALUES(?, ?)', records)
            cursor.execute("commit")
        # Create index for faster access
        cursor.execute('CREATE INDEX refseq_index ON refseq_to_uniprot (refseq_id)')
    # Remove temporal files
    if keep == False:
        if Path(input_table).is_file():
            Path(input_table).unlink()

def create_ko_to_ec(output_folder, database, keep):
    # Download information
    output_folder = Path(output_folder)
    download_output = str(output_folder / "02.KO_to_EC_identifiers.txt.gz")
    wget.download("http://enve-omics.ce.gatech.edu/data/public_microbeannotator/02.KO_to_EC_identifiers.txt.gz",
                            out=download_output)
    # Connect with database
    conn = sqlite3.connect(database)
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

    # Remove temporal files
    if keep == False:
        if Path(download_output).is_file():
            Path(download_output).unlink()
    
def create_interpro_tables(output_folder, database, keep):
    # Download information
    output_folder = Path(output_folder)
    download_output = str(output_folder / "interpro_metadata.xml.gz")
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
    conn = sqlite3.connect(database)
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
                cursor.executemany('INSERT INTO interpro_to_ec VALUES(?, ?)', records)
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
    cursor.execute('CREATE INDEX interpro_index_ce ON interpro_to_ec (interpro_id)')

    # Remove temporal files
    if keep == False:
        if Path(download_output).is_file():
            Path(download_output).unlink()

################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script build the search databases required by MicrobeAnnotator\n'''
            '''Usage: ''' + sys.argv[0] + ''' -f [output_file folder]\n'''
            '''Global mandatory parameters: -f [output_file folder]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input_table', dest='intable', action='store', required=True,
                        help='RefSeq to Uniprot table to create correspondence table.')
    parser.add_argument('-d', '--database', dest='database', action='store', required=True,
                        help='SQLite database to create tables in.')
    parser.add_argument('-o', '--output', dest='output', action='store', required=True,
                        help='Output directory where data should be stored.')
    parser.add_argument('--keep', dest='keep', action='store_true', required=False,
                        help='Keep intermediate files, can increase disk requirement (not necessary and therefore not recommended).')
    args = parser.parse_args()

    intable = args.intable
    database = args.database
    output = args.output
    keep = args.keep

    # ----------------------------
    print("Creating correspondence tables...")
    try:
        create_refseq_to_uniprot(intable, database, keep)
    except:
        print("Could not create RefSeq to Uniprot table.")
    #try:
    create_ko_to_ec(output, database, keep)
    #except:
    #    print("Could not create KO to EC table.")
    #try:
    create_interpro_tables(output, database, keep)
    #except:
    #    print("Could not create UniProt tables.")
    print("Done!")
    # ----------------------------

if __name__ == "__main__":
    main()

