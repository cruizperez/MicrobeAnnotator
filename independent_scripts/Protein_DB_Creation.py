#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      0.9
# Date:         March 05, 2020

# Description: Builds the search databases depending on the search method
# selected. Can be one of blast, sword or diamond. These must be in the PATH
# or the path must be provided.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import sys
print(sys.path)
import subprocess
from shutil import which
from sys import exit
from pathlib import Path
from .RefSeq_Data_Downloader import searching_all_files

################################################################################
"""---1.0 Define Functions---"""
def blastp_db_creator(directory, path=None):
    print('Building blast databases')
    fasta_files = searching_all_files(directory, 'protein')
    database_files = []
    if path != None:
        makeblasdb_call = Path(path) / 'makeblastdb'
    else:
        makeblasdb_call = 'makeblastdb'
    if which(makeblasdb_call) == None:
        exit(' Blast\'s "makeblastdb" not found in path, provide the folder with blast binaries')
    else:
        for protein_file in fasta_files:
            output_name = protein_file.with_suffix("")
            database_files.append(output_name)
            subprocess.call([makeblasdb_call, '-in', protein_file, '-dbtype', 'prot',
            '-parse_seqids', '-out', output_name])
    print('Done!')
    print(database_files)
    return database_files

def diamond_db_creator(directory, path=None):
    print('Building diamond databases')
    fasta_files = searching_all_files(directory, 'protein')
    database_files = []
    if path != None:
        diamond_call = Path(path) / 'makedb'
    else:
        diamond_call = 'makedb'
    if which(diamond_call) == None:
        exit('Diamond\'s "makedb" not found in path, provide the folder with diamond binaries')
    else:
        for protein_file in fasta_files:
            output_name = protein_file.with_suffix("")
            final_db_name = Path(output_name).with_suffix('.dmnd')
            database_files.append(final_db_name)
            subprocess.call([diamond_call, '--in', protein_file, '-d', output_name])
    print('Done!')
    print(database_files)
    return database_files

def sword_db_creator(directory, path=None):
    print('Building sword databases (no db needed, just checking if sword is in path)')
    fasta_files = searching_all_files(directory, 'protein')
    if path != None:
        sword_call = Path(path) / 'sword'
    else:
        sword_call = 'sword'
    if which(sword_call) == None:
        print('I did not find the sword executable. However, I don\'t need \
            it for this step, the search will probably fail in the future.\n \
            Make sure to provide the correct path for sword')
    else:
        print('Done!')
    print(fasta_files)
    return fasta_files


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script downloads the fasta and genbank files that\n'''
            '''are needed to build the RefSeq annotation database. By default it\n'''
            '''downloads both but you can specify either\n'''
            '''Usage: ''' + sys.argv[0] + ''' -f [Output folder]\n'''
            '''Global mandatory parameters: -f [Output folder]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-f', '--folder', dest='folder', action='store', required=True,
                        help='Folder where all raw fasta files are located.')
    parser.add_argument('-m', '--method', dest='method', action='store', required=True,
                        help='Search (and db creation) method, one of blast, diamond or sword')
    parser.add_argument('--bin_path', dest='bin_path', action='store', required=False,
                        help='Path to binary folder for selected method.')
    args = parser.parse_args()

    folder = args.folder
    method = args.method
    bin_path = args.bin_path
    method = method.lower()

    # ----------------------------
    if method == 'blast':
        blastp_db_creator(folder, bin_path)
    elif method == 'diamond':
        diamond_db_creator(folder, bin_path)
    elif method == 'sword':
        sword_db_creator(folder, bin_path)
    # ----------------------------

if __name__ == "__main__":
    main()
