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
from shutil import rmtree
from independent_scripts import protein_search
from independent_scripts import fasta_filter_list
from independent_scripts import sqlite3_search

################################################################################
"""---1.0 Define Functions---"""


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys, textwrap
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
                        help=textwrap.dedent('''
                        Number of processes to launch, i.e. number of protein files to process simultaneously.
                        Note this is different from threads. For more information see the README. By default 1.
                        '''))
    misc_options.add_argument('--light', dest='light', action='store_true', required=False,
                        help=textwrap.dedent('''
                        Use only KOfamscan and swissprot databases. By default also uses refseq and
                        trembl (only use if you built both using "MicrobeAnnotator_DB_Builder").
                        '''))
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
        # Final annotation fields:
        # query_id protein_id product ko_number ko_product taxonomy function compartment process database
    finally:
        pool.close()
        pool.join()
    # ----------------------------

    # Filter out the proteins already annotated
    # Create dictionary with information per protein file
    # Structure: protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it]
    print("Filtering KOFamscan results...")
    temporal_protein_folder = Path(output_dir) / "temporal_proteins"
    temporal_protein_folder.mkdir(parents=True, exist_ok=True)
    protein_file_info = {}
    for ko_result in kofam_results:
        protein_file_info[ko_result[0]] = [ko_result[1], ko_result[3]]
        outfile = str(temporal_protein_folder / (ko_result[1] + ".1it"))
        # Add first iteration protein files
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
        multiple_arguments=arguments_to_pass), input_proteins)
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

    print(protein_file_info)
    # Search annotations in SQLite DB and append to the final annotation file
    print("Extracting annotation data...")
    sql_database = str(Path(database_folder) / "02.MicrobeAnnotator.db")
    for original_file, information in protein_file_info.items():
        filtered_sword_results = information[3]
        final_annotation_file = information[1]
        significant_hits = []
        with open(filtered_sword_results) as swissprot_results:
            for line in swissprot_results:
                line = line.strip().split()
                significant_hits.append((line[0], line[1]))
        if len(significant_hits) == 0:
            continue
        else:
            annotation = sqlite3_search.search_ids_imported(sql_database, "swissprot", significant_hits)
            # Each annotation will have:
            # query_id gene_id accession product ko_number organism taxonomy function compartment process
        if light == True:
            with open(final_annotation_file, 'a') as final_annotation_fh:
                for match in annotation:
                    final_annotation_fh.write("{}\t{}\t{}\t{}\tNA\t{}\t{}\t{}\t{}\tswissprot\n".format(match[0],
                    match[1], match[3], match[4], match[6], match[7], match[8], match[9]))
            if temporal_protein_folder.is_dir():
                rmtree(temporal_protein_folder)
        else:
            # Check if all annotations have ko numbers, if not, get ids to next round of iteration
            gene_ids_no_ko = []
            with open(final_annotation_file, 'a') as final_annotation_fh:
                for match in annotation:
                    if match[4] != "":
                        final_annotation_fh.write("{}\t{}\t{}\t{}\tNA\t{}\t{}\t{}\t{}\tswissprot\n".format(match[0],
                        match[1], match[3], match[4], match[6], match[7], match[8], match[9]))
                    else:
                        gene_ids_no_ko.append(match[0])
            # If there are ids with no annotations filter the iteration 1 protein file
            if len(gene_ids_no_ko) > 0:
                second_it_outfile = str(temporal_protein_folder / (protein_file_info[original_file][0] + ".2it"))
                protein_file_info[original_file].append(second_it_outfile)
                fasta_filter_list.fastA_filter_list(protein_file_info[original_file][2], 
                second_it_outfile, gene_ids_no_ko, reverse=True)
    # -----------------------
    print(protein_file_info)

    # If running complete pipeline search against refseq
    if light == False:
        # Determine database to use
        print("Searching proteins against RefSeq...")
        if method == "blast":
            refseq_database = str(Path(database_folder) / "01.Protein_DB/refseq_protein")
        elif method == "diamond":
            refseq_database = str(Path(database_folder) / "01.Protein_DB/refseq_protein.dmnd")
        elif method == "sword":
            refseq_database = str(Path(database_folder) / "01.Protein_DB/refseq_protein.fasta")
        input_proteins = []
        for infor in protein_file_info.values():
            if len(infor) == 5:
                input_proteins.append(infor[-1])
        if len(input_proteins) > 0:
            try:
                pool = multiprocessing.Pool(processes)
                arguments_to_pass = (output_dir, 'refseq', refseq_database, method,
                                    threads, id_perc, bitscore, evalue, aln_percent, method_bin)
                search_results = pool.map(partial(protein_search.similarity_search,
                multiple_arguments=arguments_to_pass), input_proteins)
                # Results as (protein_file, filtered_search_file)
            finally:
                pool.close()
                pool.join()
        # Add name of filtered_search_results to protein_file_info
        # New structure: 
        # protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it, Swissprot_Search_File,
        # filtered_fasta_2it, RefSeq_Search_File]
        temp_dir = protein_file_info.copy()
        for result in search_results:
            for filename, info in temp_dir.items():
                if len(info) == 5 and result[0] == info[4]:
                    protein_file_info[filename].append(result[1])
        del temp_dir
        print(protein_file_info)
        # --------------------------

        # Search annotations in SQLite DB and append to the final annotation file
        for original_file, information in protein_file_info.items():
            if len(information) == 6:
                filtered_sword_results = information[5]
                final_annotation_file = information[1]
                significant_hits = []
                with open(filtered_sword_results) as swissprot_results:
                    for line in swissprot_results:
                        line = line.strip().split()
                        significant_hits.append((line[0], line[1]))
                if len(significant_hits) == 0:
                    continue
                else:
                    annotation = sqlite3_search.search_ids_imported(sql_database, "refseq", significant_hits)
                    # Each annotation will have:
                    # query_id gene_id product taxonomy ko_number ko_product
                # Check if all annotations have ko numbers, if not, get ids to next round of iteration
                gene_ids_no_ko = []
                with open(final_annotation_file, 'a') as final_annotation_fh:
                # Final annotation fields:
                # query_id protein_id product ko_number ko_product taxonomy function compartment process database
                    for match in annotation:
                        if match[4] != "" and match[4] != "NA":
                            final_annotation_fh.write("{}\t{}\t{}\t{}\t{}\t{}\tNA\tNA\tNA\trefseq\n".format(match[0],
                            match[1], match[2], match[4], match[5], match[3]))
                        else:
                            gene_ids_no_ko.append(match[0])
                # If there are ids with no annotations filter the iteration 1 protein file
                if len(gene_ids_no_ko) > 0:
                    third_it_outfile = str(temporal_protein_folder / (protein_file_info[original_file][0] + ".3it"))
                    protein_file_info[original_file].append(third_it_outfile)
                    fasta_filter_list.fastA_filter_list(protein_file_info[original_file][4], 
                    third_it_outfile, gene_ids_no_ko, reverse=True)
                    # New structure: 
                    # protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it, Swissprot_Search_File,
                    # filtered_fasta_2it, RefSeq_Search_File, filtered_fasta_3it]
        print(protein_file_info)
        # Finally run remaining proteins against Trembl
        print("Searching proteins against Trembl...")
        if method == "blast":
            refseq_database = str(Path(database_folder) / "01.Protein_DB/uniprot_trembl")
        elif method == "diamond":
            refseq_database = str(Path(database_folder) / "01.Protein_DB/uniprot_trembl.dmnd")
        elif method == "sword":
            refseq_database = str(Path(database_folder) / "01.Protein_DB/uniprot_trembl.fasta")
        input_proteins = []
        for infor in protein_file_info.values():
            if len(infor) == 7:
                input_proteins.append(infor[-1])
        if len(input_proteins) > 0:
            try:
                pool = multiprocessing.Pool(processes)
                arguments_to_pass = (output_dir, 'trembl', refseq_database, method,
                                    threads, id_perc, bitscore, evalue, aln_percent, method_bin)
                search_results = pool.map(partial(protein_search.similarity_search,
                multiple_arguments=arguments_to_pass), input_proteins)
                # Results as (protein_file, filtered_search_file)
            finally:
                pool.close()
                pool.join()
        # Add name of filtered_search_results to protein_file_info
        # New structure: 
        # protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it, Swissprot_Search_File,
        # filtered_fasta_2it, RefSeq_Search_File, filtered_fasta_3it, Trembl_Search_File]
        temp_dir = protein_file_info.copy()
        for result in search_results:
            for filename, info in temp_dir.items():
                if len(info) == 7 and result[0] == info[4]:
                    protein_file_info[filename].append(result[1])
        del temp_dir
        print(protein_file_info)
        # --------------------------

        # Search annotations in SQLite DB and append to the final annotation file
        for original_file, information in protein_file_info.items():
            if len(information) == 8:
                filtered_sword_results = information[5]
                final_annotation_file = information[1]
                significant_hits = []
                with open(filtered_sword_results) as swissprot_results:
                    for line in swissprot_results:
                        line = line.strip().split()
                        significant_hits.append((line[0], line[1]))
                if len(significant_hits) == 0:
                    continue
                else:
                    annotation = sqlite3_search.search_ids_imported(sql_database, "trembl", significant_hits)
                    # Each annotation will have:
                    # query_id gene_id product taxonomy ko_number ko_product
                # Check if all annotations have ko numbers, if not, get ids to next round of iteration
                gene_ids_no_ko = []
                with open(final_annotation_file, 'a') as final_annotation_fh:
                # Final annotation fields:
                # query_id protein_id product ko_number ko_product taxonomy function compartment process database
                    for match in annotation:
                        final_annotation_fh.write("{}\t{}\t{}\t{}\tNA\t{}\t{}\t{}\t{}\ttrembl\n".format(match[0],
                        match[1], match[3], match[4], match[6], match[7], match[8], match[9]))
        if temporal_protein_folder.is_dir():
                rmtree(temporal_protein_folder)
        print(protein_file_info)
if __name__ == "__main__":
    main()

