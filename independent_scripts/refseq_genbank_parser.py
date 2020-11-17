#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Parses compressed genbank files and extracts protein
information relevant for annotation purposes.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import argparse
from sys import argv
import gzip
import multiprocessing
import shutil
from pathlib import Path
from Bio import SeqIO

################################################################################
"""---1.0 Define Functions---"""

def table_creator(genbank_file):
    """ Parses genbank files and extract annotation information
    
    Arguments:
        genbank_file {str} -- Path to compressed genbank file
    
    Returns:
        temporal_table -- Path to temporal table created
    """
    temporal_table = Path(genbank_file).with_suffix('.temp')
    with gzip.open(genbank_file, 'rt') as uncompressed_genbank, \
    open(temporal_table, 'w') as output_file:
        for record in SeqIO.parse(uncompressed_genbank, "genbank"):
            protein_id = record.id
            product = record.description.split("[")[0].strip()
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
                    taxonomy = ", ".join(record.annotations['taxonomy']) + ", " + species
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
                if feature.type == "Protein" and "EC_number" in feature.qualifiers:
                    ec_number = feature.qualifiers["EC_number"][0]
                else:
                    ec_number = "NA"
            output_file.write("{}\t{}\t{}\t{}\n".format(protein_id, product, taxonomy, ec_number))
    return temporal_table


def table_generator_worker(genbank_list, threads):
    """ Worker calling table_creator function
    
    Arguments:
        genbank_list {list} -- List of compressed GenBank files
        threads {int} -- Threads to use
    
    Returns:
        temp_table_list -- List of temporal tables created
    """
    temp_table_list = []
    try:
        pool = multiprocessing.Pool(threads)
        temp_table_list = pool.map(table_creator, genbank_list)
    finally:
        pool.close()
        pool.join()
    return temp_table_list


def table_merger(temp_table_list, final_table, keep):
    """ Merges individual tables into a refseq annotation table
    
    Arguments:
        temp_table_list {list} -- List of paths to temporal tables created
        final_table {str} -- Final output_file table
        keep {bool} -- Retain intermediate files
    """
    with open(final_table, 'w') as output_file:
        for file in temp_table_list:
            with open(file) as temp:
                shutil.copyfileobj(temp, output_file)
            if keep == False:
                file.unlink()
            else: continue


################################################################################
"""---3.0 Main Function---"""
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script parses the downloaded compressed RefSeq Genbank files,\n'''
                        '''extracts relevant information and stores in a single table\n'''
            '''Usage: ''' + argv[0] + ''' -i [Genbank Files] -t [Threads] -o [output_file]\n'''
            '''Global mandatory parameters: -i [Genbank Files] -o [output_file]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')

    input_options = parser.add_argument_group('General input/output_file options')
    input_options.add_argument('-i', '--input_files', dest='input_files', action='store', nargs='+', required=True,
                        help='Space-separated list of compressed genbank files, or use something like "$(ls files*.gz)"')
    input_options.add_argument('-o', '--output_file', dest='output_file', action='store', required=True,
                        help='File to store final table with record information')

    misc_options = parser.add_argument_group('Miscellaneous options')
    misc_options.add_argument('-t', '--threads', dest='threads', action='store', default=1, type=int, required=False,
                        help='Number of threads to use, by default 1')
    misc_options.add_argument('-k', '--keep', dest='keep', action='store_true', required=False,
                        help='Keep intermediate files, by default false')

    args = parser.parse_args()
    
    input_files = args.input_files
    output_file = args.output_file
    threads = args.threads
    keep = args.keep

    # --------------------------------------------------
    print("Processing {} files on {} threads...".format(len(input_files), threads), end=" ")
    temp_tables = table_generator_worker(input_files, threads)
    print("Done")
    print("Merging tables...", end=" ")
    table_merger(temp_tables, output_file, keep)
    print("Done")
    # --------------------------------------------------


if __name__ == "__main__":
    main()
