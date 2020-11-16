#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Create SQLite database from a tab-delimited table.
Used to create annotation parsing databases.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import sqlite3
from sys import exit
import wget
import gzip
from pathlib import Path

################################################################################
"""---1.0 Define Functions---"""
def create_swissprot_table(database, table):
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS swissprot")
    cursor.execute('CREATE TABLE swissprot \
        (gene_id TEXT, accession TEXT, gene_name TEXT, ko_number TEXT, \
        organism TEXT, taxonomy TEXT, function TEXT, compartment TEXT, \
        process TEXT, interpro TEXT, pfam TEXT, ec_number TEXT)')
    with open(table, 'r') as swissprot:
        record_counter = 0
        records = []
        for line in swissprot:
            # Commit changes after 500000 records
            if record_counter == 500000:
                cursor.execute("begin")
                cursor.executemany('INSERT INTO swissprot VALUES(?, ?, ?, ?, ?, \
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
            cursor.executemany('INSERT INTO swissprot VALUES(?, ?, ?, ?, ?, \
                ?, ?, ?, ?, ?, ?, ?)', records)
            cursor.execute("commit")
        # Create index for faster access
        cursor.execute('CREATE INDEX swissprot_id ON swissprot (gene_id)')

def create_trembl_table(database, table):
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS trembl")
    cursor.execute('CREATE TABLE trembl \
        (gene_id TEXT, accession TEXT, gene_name TEXT, ko_number TEXT, \
        organism TEXT, taxonomy TEXT, function TEXT, compartment TEXT, \
        process TEXT, interpro TEXT, pfam TEXT, ec_number)')
    with open(table, 'r') as trembl:
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

def create_refseq_table(database, table):
    # Store protein ko annotations
    # Download refseq KO annotations
    print("Parsing RefSeq KO numbers")
    annotation_file = wget.download("http://enve-omics.ce.gatech.edu/data/public_microbeannotator/01.RefSeq_KEGG_Annotations.txt.gz",
                            out="refseq_ko_annotations.tar.gz")
    protein_annotations = {}
    with gzip.open(annotation_file, 'rt') as annotations:
        for line in annotations:
            line = line.rstrip('\n').split("\t")
            protein_annotations[line[0]] = [line[1], line[2]]
    Path(annotation_file).unlink()
    print("Done")
    # Connect to db and create table
    print("Creating RefSeq table")
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS refseq")
    cursor.execute('CREATE TABLE refseq \
        (gene_id TEXT, product TEXT, taxonomy TEXT, \
        ec_number TEXT, ko_number TEXT, ko_product TEXT)')
    with open(table, 'r') as refseq:
        record_counter = 0
        records = []
        for line in refseq:
            # Commit changes after 500000 records
            if record_counter == 500000:
                cursor.execute("begin")
                cursor.executemany('INSERT INTO refseq VALUES(?, ?, ?, ?, ?, ?)', records)
                cursor.execute("commit")
                record_counter = 0
                records = []
            else:
                line = line.rstrip('\n').split("\t")
                if line[0] in protein_annotations:
                    ko_number = protein_annotations[line[0]][0]
                    ko_product = protein_annotations[line[0]][1]
                    records.append((line[0], line[1], line[2], line[3], ko_number, ko_product))
                    record_counter += 1
                else:
                    records.append((line[0], line[1], line[2], line[3], "NA", "NA"))
                    record_counter += 1
        # Commit remaining records
        if record_counter > 0:
            cursor.execute("begin")
            cursor.executemany('INSERT INTO refseq VALUES(?, ?, ?, ?, ?, ?)', records)
            cursor.execute("commit")
        # Create index for faster access
        cursor.execute('CREATE INDEX refseq_id ON refseq (gene_id)')
    print("Done")


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script builds a sqlite database from a tab separated table.\n'''
            '''By default it assumes the first line of the input table has headers, if not\n'''
            '''you must provide a list of headers as --header header1,header2,header3...\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [Input Table] -d [Database Name]\n'''
            '''Global mandatory parameters: -i [Input Table] -d [Database Name]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_table', action='store', required=True,
                        help='Input tab-delimited table to parse, by default assumes headers are present')
    parser.add_argument('-d', '--database', dest='database', action='store', required=True,
                        help='Database name, can be exisiting or new.')
    parser.add_argument('--type', dest='db_type', action='store', required=True,
                        help='Table type, can be one of swissprot, trembl, or refseq')
    args = parser.parse_args()

    input_table = args.input_table
    database = args.database
    db_type = args.db_type
    db_type = db_type.lower()

    # ----------------------------
    print("Building database")
    if db_type == "swissprot":
        create_swissprot_table(database, input_table)
    elif db_type == "trembl":
        create_trembl_table(database, input_table)
    elif db_type == "refseq":
        create_refseq_table(database, input_table)
    else:
        exit("Please specify a table type, swissprot, trembl or refseq")
    print("Database succesfully created!")
    # ----------------------------

if __name__ == "__main__":
    main()
