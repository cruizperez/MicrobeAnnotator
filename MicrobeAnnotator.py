#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      0.9
# Date:         May 06, 2020

# Description: Downloads all data required to build the search databases.
# Parses the annotation information and creates the method-specific dbs.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import multiprocessing
from functools import partial
from pathlib import Path
from independent_scripts import protein_search
from independent_scripts import fasta_filter_list

################################################################################
"""---1.0 Define Functions---"""


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
    general_options = parser.add_argument_group('Mandatory i/o options.')
    general_options.add_argument('-i', '--input', dest='input_list', action='store', required=True, nargs='+',
                        help='Space-separated list of protein files to parse.')
    general_options.add_argument('-o', '--outdir', dest='output_dir', action='store', required=True,
                        help='Directory to store results.')
    general_options.add_argument('-d', '--database', dest='database_folder', action='store', required=True,
                        help='Directory where MicrobeAnnotator databases are located.')
    search_options = parser.add_argument_group('Options for search process.')
    search_options.add_argument('-m', '--method', dest='method', action='store', required=True,
                        help='Method used to create databases and to perform seaches. One of "blast", "diamond" or "sword".')
    search_options.add_argument('--kofam_bin', dest='kofam_bin', action='store', required=False, default=None,
                        help='Directory where KOFamscan binaries are located. By default assumes it is in PATH.')
    search_options.add_argument('--method_bin', dest='method_bin', action='store', required=False, default=None,
                        help='Directory where KOFamscan binaries are located. By default assumes it is in PATH.')
    search_options.add_argument('--id_perc', dest='id_perc', action='store', required=False, default=40, type=int,
                        help='Minimum identity percentage to retain a hit. By default 40.')
    search_options.add_argument('--bitscore', dest='bitscore', action='store', required=False, default=50, type=int,
                        help='Minimum bitscore to retain a hit. By default 50.')
    search_options.add_argument('--evalue', dest='evalue', action='store', required=False, default=0.01, type=float,
                        help='Maximum evalue to retain a hit. By default 0.01.')
    search_options.add_argument('--aln_percent', dest='aln_percent', action='store', required=False, default=60, type=int,
                        help='Minimum percentage of query covered by hit alignment. By default 60.')
    misc_options = parser.add_argument_group('Miscellaneous options')
    misc_options.add_argument('-t', '--threads', dest='threads', action='store', required=False, default=1, type=int,
                        help='Threads to use per processed file, i.e. (per protein file). By default 1.')
    misc_options.add_argument('-p', '--processes', dest='processes', action='store', required=False, default=1, type=int,
                        help='Number of processes to launch, i.e. number of protein files to process simultaneously. \
                            Note this is different from threads. For more information see the README. By default 1.')
    misc_options.add_argument('--light', dest='light', action='store_true', required=False,
                        help='Use only KOfamscan and swissprot databases. By default also uses refseq and trembl (only use if you built both using "MicrobeAnnotator_DB_Builder").')
    args = parser.parse_args()

    input_list = args.input_list
    output_dir = args.output_dir
    database_folder = args.database_folder
    method = args.method
    kofam_bin = args.kofam_bin
    method_bin = args.method_bin
    id_perc = args.id_perc
    bitscore = args.bitscore
    evalue = args.evalue
    aln_percent = args.aln_percent
    threads = args.threads
    processes = args.processes
    light = args.light


    # ----------------------------
    # Search the initial dataset with KOFamscan
    print("Searching proteins using KOFamscan...")
    try:
        pool = multiprocessing.Pool(processes)
        kofam_results = pool.map(partial(protein_search.kofamscan_annotation, 
        multi_argument=(output_dir, threads, kofam_bin)), input_list)
        # Results as (protein_file, protein_file_name, ids_proteins_annotated, final_file)
    finally:
        pool.close()
        pool.join()
    # ----------------------------

    # Filter out the proteins already annotated
    # Create dictionary with information per protein file
    # Structure: protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it]
    print("Filtering KOFamscan results...")
    temporal_protein_folder = Path(output_directory) / "temporal_proteins"
    temporal_protein_folder.mkdir(parents=True, exist_ok=True)
    protein_file_info = {}
    for ko_result in kofam_results:
        protein_file_info[ko_result[0]] = [ko_result[1], ko_result[3]]
        outfile = str(temporal_protein_folder / ko_result[1] + ".1it.faa")
        protein_file_info[ko_result[0]].append(outfile)
        fasta_filter_list.fastA_filter_list(ko_result[0], outfile, ko_result[2], reverse=True)
    # ----------------------------

    # Search proteins NOT annotated with KOFamscan against Swissprot
    # Determine database to use
    print("Searching proteins against Swissprot...")
    if method == "blast":
        swissprot_database = str(Path(database_folder) / "01.Protein_DB/uniprot_sprot")
    elif method == "diamond":
        swissprot_database = str(Path(database_folder) / "01.Protein_DB/uniprot_sprot.dmnd")
    elif method == "sword":
        swissprot_database = str(Path(database_folder) / "01.Protein_DB/uniprot_sprot.fasta")
    input_proteins = []
    for infor in protein_file_info.values():
        input_proteins.append(infor[-1])
    try:
        pool = multiprocessing.Pool(processes)
        arguments_to_pass = (output_dir, 'swissprot', swissprot_database, method,
                            threads, id_perc, bitscore, evalue, aln_percent, method_bin)
        search_results = pool.map(partial(protein_search.similarity_search,
        multi_argument=arguments_to_pass), input_proteins)
        # Results as (protein_file, filtered_search_file)
    finally:
        pool.close()
        pool.join()
    # Add name of filtered_search_results to protein_file_info
    # New structure: 
    # protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it, Swissprot_Search_File]
    temp_dir = protein_file_info.copy()
    for result in search_results:
        for filename, info in temp_dir.items():
            if result[0] == info[2]:
                protein_file_info[filename].append(result[1])
    del temp_dir
    # ----------------------------
    #TODO buscar en sql db hits y añadir resultados a annotation file (sqlite_search - END)
    #TODO filtrar aquellos que no tengan KO (If light == false)
    #TODO repetir para RefSeq y Trembl
    #TODO sacar todos los KO per file y mappearlos (ko_mapper)

if __name__ == "__main__":
    main()

