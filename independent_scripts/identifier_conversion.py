#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Finds identifiers from one database into another.
# Uses the interconversion tables in the main database.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import sqlite3

################################################################################
"""---1.0 Define Functions---"""
def convert_refseq_to_uniprot(identifier_list, sql_database, inverse):
    conn = sqlite3.connect(sql_database)
    cursor = conn.cursor()
    refseq_to_uniprot = {}
    uniprot_to_refseq = {}
    for identifier in identifier_list:
        # If uniprot to refseq
        if inverse == False:
            cursor.execute("SELECT * FROM refseq_to_uniprot WHERE uniprot_id=?", (identifier,))
            rows = cursor.fetchall()
            if len(rows) > 0:
                uniprot_to_refseq[identifier] = []
                for match in rows:
                    # Append refseq matches
                    uniprot_to_refseq[identifier].append(match[0])
            return uniprot_to_refseq
        # Else search refseq to uniprot
        else:
            cursor.execute("SELECT * FROM refseq_to_uniprot WHERE refseq_id=?", (identifier,))
            rows = cursor.fetchall()
            if len(rows) > 0:
                refseq_to_uniprot[identifier] = []
                for match in rows:
                    # Append uniprot matches
                    refseq_to_uniprot[identifier].append(match[1])
            return refseq_to_uniprot


def convert_ko_to_ec(identifier_list, sql_database, inverse):
    conn = sqlite3.connect(sql_database)
    cursor = conn.cursor()
    ko_to_ec = {}
    ec_to_ko = {}
    for identifier in identifier_list:
        # If ec to ko
        if inverse == False:
            cursor.execute("SELECT * FROM ko_to_ec WHERE ec_identifier=?", (identifier,))
            rows = cursor.fetchall()
            if len(rows) > 0:
                ec_to_ko[identifier] = []
                for match in rows:
                    # Append ko matches
                    ec_to_ko[identifier].append(match[0])
            return ec_to_ko
        # Else search ko to ec
        else:
            cursor.execute("SELECT * FROM ko_to_ec WHERE ko_identifier=?", (identifier,))
            rows = cursor.fetchall()
            if len(rows) > 0:
                ko_to_ec[identifier] = []
                for match in rows:
                    # Append ec matches
                    ko_to_ec[identifier].append(match[1])
            return ko_to_ec

def convert_interpro_to_ko(identifier_list, sql_database, inverse):
    conn = sqlite3.connect(sql_database)
    cursor = conn.cursor()
    interpro_to_ko = {}
    ko_to_interpro = {}
    for identifier in identifier_list:
        # If ko to interpro
        if inverse == True:
            cursor.execute("SELECT * FROM interpro_to_ko WHERE ko_identifier=?", (identifier,))
            rows = cursor.fetchall()
            if len(rows) > 0:
                ko_to_interpro[identifier] = []
                for match in rows:
                    # Append ko matches
                    ko_to_interpro[identifier].append(match[0])
            return ko_to_interpro
        # Else search interpro to ko
        else:
            cursor.execute("SELECT * FROM interpro_to_ko WHERE interpro_id=?", (identifier,))
            rows = cursor.fetchall()
            if len(rows) > 0:
                interpro_to_ko[identifier] = []
                for match in rows:
                    # Append ec matches
                    interpro_to_ko[identifier].append(match[1])
            return interpro_to_ko

def convert_interpro_to_ec(identifier_list, sql_database, inverse):
    conn = sqlite3.connect(sql_database)
    cursor = conn.cursor()
    interpro_to_ec = {}
    ec_to_interpro = {}
    for identifier in identifier_list:
        # If ko to interpro
        if inverse == True:
            cursor.execute("SELECT * FROM interpro_to_ec WHERE ec_identifier=?", (identifier,))
            rows = cursor.fetchall()
            if len(rows) > 0:
                ec_to_interpro[identifier] = []
                for match in rows:
                    # Append ko matches
                    ec_to_interpro[identifier].append(match[0])
            return ec_to_interpro
        # Else search interpro to ko
        else:
            cursor.execute("SELECT * FROM interpro_to_ec WHERE interpro_id=?", (identifier,))
            rows = cursor.fetchall()
            if len(rows) > 0:
                interpro_to_ec[identifier] = []
                for match in rows:
                    # Append ec matches
                    interpro_to_ec[identifier].append(match[1])
            return interpro_to_ec


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script searches for identifiers and returns associated ids from other databases.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -f [output_file folder]\n'''
            '''Global mandatory parameters: -f [output_file folder]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    general_options = parser.add_argument_group('Inputs and outputs.')
    general_options.add_argument('-i', '--input_file', dest='input_file', action='store', required=True,
                        help='File with identifiers to search and convert (one per line).')
    general_options.add_argument('-d', '--database', dest='database', action='store', required=True,
                        help='SQLite database to search.')
    general_options.add_argument('-o', '--outfile', dest='outfile', action='store', required=True,
                        help='Output table with original id and corresponding ids in other databases.')
    conversion_to_perform = parser.add_argument_group('Available conversions.')
    conversion_to_perform.add_argument('--refseq_to_uniprot', dest='refseq_to_uniprot', action='store_true', required=False,
                        help='Convert refseq -> uniprot or uniprot -> refseq (see --inverse)')
    conversion_to_perform.add_argument('--ko_to_ec', dest='ko_to_ec', action='store_true', required=False,
                        help='Convert ko -> ec or ec -> ko (see --inverse)')
    conversion_to_perform.add_argument('--interpro_to_ko', dest='interpro_to_ko', action='store_true', required=False,
                        help='Convert interpro -> ko or ko -> interpro (see --inverse)')
    conversion_to_perform.add_argument('--interpro_to_ec', dest='interpro_to_ec', action='store_true', required=False,
                        help='Convert interpro -> ec or ec -> interpro (see --inverse)')
    conversion_modificatopn = parser.add_argument_group('Changes to searching order.')
    conversion_modificatopn.add_argument('--inverse', dest='inverse', action='store_true', required=False, default=False,
                        help='Change the order of conversion. For example, in ko_to_ec instead to converting ko -> ec, convert ec -> ko.')
    args = parser.parse_args()

    input_file = args.input_file
    database = args.database
    outfile = args.outfile
    refseq_to_uniprot = args.refseq_to_uniprot
    ko_to_ec = args.ko_to_ec
    interpro_to_ko = args.interpro_to_ko
    interpro_to_ec = args.interpro_to_ec
    inverse = args.inverse

    # ----------------------------
    print("Checking inputs...")
    # Check if more than one or no conversion was selected
    conversion_types = [refseq_to_uniprot, ko_to_ec, interpro_to_ko, interpro_to_ec]
    if conversion_types.count(True) == 0:
        exit("Plase select one type of conversion to perform.")
    elif conversion_types.count(True) > 1:
        exit("Plase select only one type of conversion to perform.")
    
    # Parse the input file and get the identifiers
    identifier_list = []
    with open(input_file, 'r') as infile:
        for line in infile:
            identifier_list.append(line.strip().split()[0])
    print("Done!")

    print("Searching in the database...")
    corresponding_ids = None
    if refseq_to_uniprot == True:
        corresponding_ids = convert_refseq_to_uniprot(identifier_list, database, inverse)
    elif ko_to_ec == True:
        corresponding_ids = convert_ko_to_ec(identifier_list, database, inverse)
    elif interpro_to_ko == True:
        corresponding_ids = convert_interpro_to_ko(identifier_list, database, inverse)
    elif interpro_to_ec == True:
        corresponding_ids = convert_interpro_to_ec(identifier_list, database, inverse)
    print("Done!")

    print("Writing output...")
    with open(outfile, 'w') as output:
        for identifier, matches in corresponding_ids.items():
            for match in matches:
                output.write("{}\t{}\n".format(identifier, match))
    print("Done!")
    # ----------------------------

if __name__ == "__main__":
    main()