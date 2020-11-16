#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Builds the search databases depending on the search method
# selected. Can be one of blast, sword or diamond. These must be in the PATH
# or the path must be provided.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import subprocess
from shutil import which
from sys import exit
from pathlib import Path

################################################################################
"""---1.0 Define Functions---"""
def blastp_db_creator(directory, path=None):
    print('Building blast databases')
    database_files = []
    if path != None:
        makeblasdb_call = Path(path) / 'makeblastdb'
    else:
        makeblasdb_call = 'makeblastdb'
    if which(makeblasdb_call) == None:
        script_path = Path(__file__)
        script_dir = script_path.parent
        stand_alone_script = script_dir / "protein_db_creation.py"
        print('Blast\'s "makeblastdb" not found in path, I cannot continue and so, I will exit.')
        print('If you are running the microbeannotator_db_builder script this is the last step of the process.')
        print('Given that this last step failed you need to run the standalone ' + str(stand_alone_script) + ' script')
        print('and provide the binary folder for the selected method using the "--bin_path" option or add')
        print('the binary to your PATH (see README for help on this).')
        exit('Otherwise, just re-run this script providing the appropriate path. :)')
    else:
        for protein_file in Path(directory).iterdir():
            if protein_file.suffix == '.fasta':
                output_file_name = protein_file.with_suffix("")
                database_files.append(output_file_name)
                subprocess.call([makeblasdb_call, '-in', protein_file, '-dbtype', 'prot',
                '-out', output_file_name])
    print('Done!')
    with open(Path(directory)/'blast_db.list', 'w') as db_list:
        for file in database_files:
            db_list.write('{}\n'.format(file.name))

def diamond_db_creator(directory, threads, path=None):
    print('Building diamond databases')
    database_files = []
    if path != None:
        diamond_call = Path(path) / 'diamond'
    else:
        diamond_call = 'diamond'
    if which(diamond_call) == None:
        script_path = Path(__file__)
        script_dir = script_path.parent
        stand_alone_script = script_dir / "protein_db_creation.py"
        print('Diamond\'s "makedb" not found in path, I cannot continue and so, I will exit.')
        print('If you are running the microbeannotator_db_builder script this is the last step of the process.')
        print('Given that this last step failed you need to run the standalone ' + str(stand_alone_script) + ' script')
        print('and provide the binary folder for the selected method using the "--bin_path" option or add')
        print('the binary to your PATH (see README for help on this).')
        exit('Otherwise, just re-run this script providing the appropriate path. :)')
    else:
        for protein_file in Path(directory).iterdir():
            if protein_file.suffix == '.fasta':
                output_file_name = protein_file.with_suffix("")
                final_db_name = Path(output_file_name).with_suffix('.dmnd')
                database_files.append(final_db_name)
                subprocess.call([diamond_call, 'makedb', '--in', protein_file, '-d', output_file_name, '--threads', str(threads)])
    print('Done!')
    with open(Path(directory)/'diamond_db.list', 'w') as db_list:
        for file in database_files:
            db_list.write('{}\n'.format(file.name))

def sword_db_creator(directory, path=None):
    print('Building sword databases (no db needed, just checking if sword is in path)')
    if path != None:
        sword_call = Path(path) / 'sword'
    else:
        sword_call = 'sword'
    if which(sword_call) == None:
        print('I did not find the sword executable. I don\'t need')
        print('it for this step. However, the search might fail in the future.')
        print('Make sure to provide the correct path for sword when searching.')
    else:
        print('Done!')
    with open(Path(directory)/'sword_db.list', 'w') as db_list:
        for file in Path(directory).iterdir():
            if file.suffix == '.fasta':
                db_list.write('{}\n'.format(file.name))


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script downloads the fasta and genbank files that\n'''
            '''are needed to build the RefSeq annotation database. By default it\n'''
            '''downloads both but you can specify either\n'''
            '''Usage: ''' + sys.argv[0] + ''' -f [output_file folder]\n'''
            '''Global mandatory parameters: -f [output_file folder]\n'''
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
