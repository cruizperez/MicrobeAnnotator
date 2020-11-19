#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: MicrobeAnnotator parses protein fasta files and annotates them
# using several databases in an iterative fashion and summarizes the findings
# using KEGG modules based on KO numbers associated with best database matches.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
from independent_scripts import identifier_conversion
import multiprocessing
from functools import partial
from pathlib import Path
from shutil import rmtree
from shutil import copyfile
from independent_scripts import protein_search
from independent_scripts import fasta_filter_list
from independent_scripts import sqlite3_search
from independent_scripts import ko_mapper
from shutil import which
import pickle

################################################################################
"""---1.0 Main Function---"""

def main():
    import argparse, sys, textwrap
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''MicrobeAnnotator parses protein fasta files and annotates them\n'''
            '''using several databases in an iterative fashion and summarizes the findings\n'''
            '''using KEGG modules based on KO numbers associated with best database matches.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [protein file] -o [output folder] -d [MicrobeAnnotator db folder]\n'''
            '''-m [search method]\n'''
            '''Global mandatory parameters: -i [protein file] -o [output folder] -d [MicrobeAnnotator db folder]\n'''
            '''-m [search method]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    general_options = parser.add_argument_group('Mandatory i/o options.')
    general_options.add_argument('-i', '--input', dest='input_list', action='store', required=False, nargs='+',
                        help='Space-separated list of protein files to parse. Use -i OR -l.')
    general_options.add_argument('-l', '--list', dest='file_list', action='store', required=False,
                        help='File with list of inputs. Use -i OR -l.')
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
    search_options.add_argument('--bitscore', dest='bitscore', action='store', required=False, default=80, type=int,
                        help='Minimum bitscore to retain a hit. By default 80.')
    search_options.add_argument('--evalue', dest='evalue', action='store', required=False, default=0.01, type=float,
                        help='Maximum evalue to retain a hit. By default 0.01.')
    search_options.add_argument('--aln_percent', dest='aln_percent', action='store', required=False, default=70, type=int,
                        help='Minimum percentage of query covered by hit alignment. By default 70.')
    plot_options = parser.add_argument_group('Summary and plotting options.')
    plot_options.add_argument('--cluster', dest='cluster', action='store', required=False,
                        help=textwrap.dedent('''
                        Cluster genomes and/or modules. Select "cols" for genomes, "rows" for modules, or "both".
                        By default, no clustering
                        '''))
    plot_options.add_argument('--filename', dest='plot_filename', action='store', required=False, default='metabolic_summary_',
                        help='Prefix for output summary tables and plots. By default "metabolic_summary"')
    misc_options = parser.add_argument_group('Miscellaneous options.')
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
                        trembl (only use if you built both using "microbeannotator_db_builder").
                        '''))
    misc_options.add_argument('--full', dest='full', action='store_true', required=False,
                        help=textwrap.dedent('''
                        Do not perform the iterative annotation but search all proteins against all databases
                        (Increases computation time).
                        '''))
    misc_options.add_argument('--continue_run', dest='continue_run', action='store_true', required=False,
                        help=textwrap.dedent('''
                        If something went wrong when runnin MicrobeAnnotator, try to resume from the last completed step.
                        '''))
    misc_options.add_argument('--refine', dest='refine', action='store_true', required=False,
                        help=textwrap.dedent('''
                        Complement the annotations by finding links to identifiers in other databases.
                        '''))
    args = parser.parse_args()

    input_list = args.input_list
    file_list = args.file_list
    output_dir = args.output_dir
    database_folder = args.database_folder
    method = args.method
    method = method.lower()
    kofam_bin = args.kofam_bin
    method_bin = args.method_bin
    id_perc = args.id_perc
    bitscore = args.bitscore
    evalue = args.evalue
    aln_percent = args.aln_percent
    threads = args.threads
    processes = args.processes
    light = args.light
    cluster = args.cluster
    if cluster != None:
        cluster = cluster.lower()
    plot_filename = args.plot_filename
    full = args.full
    refine = args.refine
    continue_run = args.continue_run

    # Check user input
    if method != 'blast' and method != 'diamond' and method != 'sword':
        exit('Please provide a valid search method with "-m", one of "blast", "diamond" or "sword"')
    print("---- This is MicrobeAnnotator ----")
    if input_list != None:
        print('I will process a total of {} files. Of these, I will run {} in parallel with {} threads per file.'.format(
            str(len(input_list)), str(processes), str(threads)))
    elif file_list != None:
        number_files = 0
        with open(file_list) as file_list_handler:
            for line in file_list_handler:
                number_files += 1
        print('I will process a total of {} files. Of these, I will run {} in parallel with {} threads per file.'.format(
            str(number_files), str(processes), str(threads)))
    if light == True:
        print('For this run I will use KOFamscan and the Swissprot database to annotate your proteins.')
    else:
        print('For this run I will use KOFamscan and the Swissprot, RefSeq and Trembl databases to annotate your proteins.')
    if cluster == None:
        print('When I\'m done I will plot the metabolism summaries using {} as file name prefix using no clustering.'.format(plot_filename))
    elif cluster == 'cols':
        print('When I\'m done I will plot the metabolism summaries using {} as file name prefix and cluster the genomes (cols)'.format(plot_filename))
    elif cluster == 'rows':
        print('When I\'m done I will plot the metabolism summaries using {} as file name prefix and cluster the modules (rows)'.format(plot_filename))
    elif cluster == 'both':
        print('When I\'m done I will plot the metabolism summaries using {} as file name prefix and cluster the genomes and modules'.format(plot_filename))
    else:
        exit('I do not recognize this clustering option: {}.'.format(cluster))
    print("---------\n")

    # Check kofamscan binary existence
    if kofam_bin == None:
        kofam_call = 'exec_annotation'
    else:
        kofam_call = str(Path(kofam_bin) / 'exec_annotation')
    if which(kofam_call) == None:
            if kofam_bin == None:
                exit("KOfamscan not found in PATH, please provide the correct path to the folder with binaries")
            else:
                exit("KOfamscan not found in " + kofam_bin + ", please provide the correct path to the folder with binaries")
    
    # Check search method binary existence
    call_name = None
    if method == 'blast':
        call_name = 'blastp'
    elif method == 'diamond':
        call_name = 'diamond'
    elif method == 'sword':
        call_name = 'sword'
    if method_bin == None:
        method_call = call_name
    else:
        method_call = str(Path(method_bin) / call_name)
    if which(method_call) == None:
            if method_bin == None:
                exit(method_call + " not found in PATH, please provide the correct path to the folder with binaries")
            else:
                exit(method_call + " not found in " + method_bin + ", please provide the correct path to the folder with binaries")
    # ----------------------------

    # Get the initial list of proteins per genome and store them
    if input_list != None and file_list != None:
        exit('Please provide only a file with the list of inputs "-l" OR a space-separated list of files "-i".')
    elif file_list != None:
        input_list = []
        with open(file_list) as file_list_handler:
            for line in file_list_handler:
                input_list.append(line.strip())
    starting_proteins = {}
    for file in input_list:
        starting_filename = str(Path(file).name)
        starting_proteins[starting_filename] = []
        with open(Path(file)) as proteins:
            for line in proteins:
                if line.startswith(">"):
                    line = line.strip().split()[0].replace(">", "")
                    starting_proteins[starting_filename].append(line)
    # ----------------------------

    # Create log folder to track process and folder to store temporal proteins
    process_step = 0
    process_log_folder = Path(output_dir) / "process_log"
    process_log_folder.mkdir(parents=True, exist_ok=True)
    temporal_protein_folder = Path(output_dir) / "temporal_proteins"
    temporal_protein_folder.mkdir(parents=True, exist_ok=True)
    sql_database = str(Path(database_folder) / "02.MicrobeAnnotator.db")
    interconversion_database = str(Path(database_folder) / "03.Conversion.db")
    # Create dictionary with information per protein file
    protein_file_info = {}
    # ----------------------------

    # Check if continue_run was passed, and if so, import the last step completed and the dictionary with information
    if continue_run == True:
        print("Restarting from last checkpoint")
        with open(process_log_folder / "log.txt", 'r') as logfile:
            process_step = int(logfile.readline().strip())
            print(process_step)
        with open(process_log_folder / "structure.pickle", 'rb') as structure_file:
            protein_file_info = pickle.load(structure_file)
    
    if process_step == 0:
        # Search the initial dataset with KOFamscan
        print("Searching proteins using KOFamscan...")
        try:
            pool = multiprocessing.Pool(processes)
            kofam_results = pool.map(partial(protein_search.kofamscan_annotation, 
            multi_argument=(output_dir, threads, kofam_bin)), input_list)
            # Results as (protein_file, protein_file_name, ids_proteins_annotated, final_file)
            # Final annotation fields:
            # * query_id protein_id product ko_number ko_product taxonomy function_go compartment_go process_go interpro pfam ec_number database
        finally:
            pool.close()
            pool.join()
        # ----------------------------

        # Filter out the proteins already annotated
        
        print("Filtering KOFamscan results...")
        # New structure after first step: protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it]
        for ko_result in kofam_results:
            protein_file_info[ko_result[0]] = [ko_result[1], ko_result[3]]
            outfile = str(temporal_protein_folder / (ko_result[1] + ".1it"))
            # Add first iteration protein files
            protein_file_info[ko_result[0]].append(outfile)
            if full == True:
                copyfile(ko_result[0],outfile)
            else:
                for protein in ko_result[2]:
                    starting_proteins[str(Path(ko_result[0]).name)].remove(protein)
                fasta_filter_list.fastA_filter_list(ko_result[0], outfile, ko_result[2], reverse=True)
        # When the process is complete, write step completed (+1) in the log folder
        # Also, export dictionary with information to be imported in case of continue
        process_step += 1
        with open(process_log_folder / "log.txt", 'w') as logfile:
            logfile.write("{}".format(process_step))
        with open(process_log_folder / "structure.pickle", 'wb') as structure_file:
            pickle.dump(protein_file_info, structure_file)
        # ----------------------------

    if process_step == 1:
        # Search proteins NOT annotated with KOFamscan against Swissprot
        # Determine database to use
        print("Searching proteins against Swissprot using " + method + " ...")
        swissprot_database = None
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

        # Search annotations in SQLite DB and append to the final annotation file
        print("Extracting Swissprot annotation data...")
        for original_file, information in protein_file_info.items():
            filtered_sword_results = information[3]
            final_annotation_file = information[1]
            significant_hits = []
            with open(filtered_sword_results) as swissprot_results:
                for line in swissprot_results:
                    line = line.strip().split()
                    significant_hits.append((line[0], line[1]))
            if len(significant_hits) == 0:
                # If no hits were found
                if light == True:
                    # Write all the proteins without annotation
                    with open(final_annotation_file, 'a') as final_annotation_fh:
                        for original_protein in starting_proteins[str(Path(original_file).name)]:
                            final_annotation_fh.write("{}\tNA\tNo match found\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n".format(original_protein))
                else:
                    # Copy the first iteration proteins to the second iteration file
                    second_it_outfile = str(temporal_protein_folder / (protein_file_info[original_file][0] + ".2it"))
                    copyfile(protein_file_info[original_file][2],second_it_outfile)
                    protein_file_info[original_file].append(second_it_outfile)
            else:
                # Extract the annotations
                annotation = sqlite3_search.search_ids_imported(sql_database, "swissprot", significant_hits)
                # Each annotation will have:
                # query_id gene_id accession product ko_number organism taxonomy function_GO compartment_GO process_GO interpro_id pfam_id EC_number
                if light == True:
                    # Write the annotations found in the annotation file
                    with open(final_annotation_file, 'a') as final_annotation_fh:
                        for match in annotation:
                            final_annotation_fh.write("{}\t{}\t{}\t{}\tNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tswissprot\n".format(match[0],
                            match[1], match[3], match[4], match[6], match[7], match[8], match[9], match[10], match[11], match[12]))
                            if match[4] != "":
                                starting_proteins[str(Path(original_file).name)].remove(match[0])
                    # Check which proteins were not annotated and add information on those
                    # Extract annotated proteins
                    with open(final_annotation_file, 'a') as final_annotation_fh:
                        for original_protein in starting_proteins[str(Path(original_file).name)]:
                            final_annotation_fh.write("{}\tNA\tNo match found\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n".format(original_protein))
                else:
                    # Write the annotations found in the annotation file
                    with open(final_annotation_file, 'a') as final_annotation_fh:
                        for match in annotation:
                            final_annotation_fh.write("{}\t{}\t{}\t{}\tNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tswissprot\n".format(match[0],
                            match[1], match[3], match[4], match[6], match[7], match[8], match[9], match[10], match[11], match[12]))
                            if match[4] != "" and match[4] != "NA":
                                starting_proteins[str(Path(original_file).name)].remove(match[0])
                    # Check which proteins were not annotated and filter for the next iteration
                    second_it_outfile = str(temporal_protein_folder / (protein_file_info[original_file][0] + ".2it"))
                    if full == True:
                        copyfile(protein_file_info[original_file][2],second_it_outfile)
                        protein_file_info[original_file].append(second_it_outfile)
                    else:
                        protein_file_info[original_file].append(second_it_outfile)
                        to_retain = starting_proteins[str(Path(original_file).name)].copy()
                        fasta_filter_list.fastA_filter_list(protein_file_info[original_file][2], 
                        second_it_outfile, to_retain, reverse=False)
                        # starting_proteins[str(Path(original_file).name)]
                    # New structure: 
                    # protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it, Swissprot_Search_File,
                    # filtered_fasta_2it]
        # When the process is complete, write step completed (+1) in the log folder
        # Also, export dictionary with information to be imported in case of continue
        process_step += 1
        with open(process_log_folder / "log.txt", 'w') as logfile:
            logfile.write("{}".format(process_step))
        with open(process_log_folder / "structure.pickle", 'wb') as structure_file:
            pickle.dump(protein_file_info, structure_file)
        # -----------------------
    
    
    
    # If running complete pipeline search against refseq
    if process_step == 2:
        if light == True:
            print("The light run of MicrobeAnnotator had finished succesfully!")
            print("Summarizing those results...")
            process_step = 5
        else:
            # Determine database to use
            print("Searching proteins against RefSeq...")
            refseq_database = None
            if method == "blast":
                refseq_database = str(Path(database_folder) / "01.Protein_DB/refseq_protein")
            elif method == "diamond":
                refseq_database = str(Path(database_folder) / "01.Protein_DB/refseq_protein.dmnd")
            elif method == "sword":
                refseq_database = str(Path(database_folder) / "01.Protein_DB/refseq_protein.fasta")
            # If we have a second iteration file then perform the searches againts refseq
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
            # --------------------------

            # Search annotations in SQLite DB and append to the final annotation file
            print("Extracting RefSeq annotation data...\n")
            for original_file, information in protein_file_info.items():
                if len(information) == 6:
                    filtered_sword_results = information[5]
                    final_annotation_file = information[1]
                    significant_hits = []
                    with open(filtered_sword_results) as swissprot_results:
                        for line in swissprot_results:
                            line = line.strip().split()
                            significant_hits.append((line[0], line[1]))
                    # If no hits were found
                    if len(significant_hits) == 0:
                        third_it_outfile = str(temporal_protein_folder / (protein_file_info[original_file][0] + ".3it"))
                        copyfile(protein_file_info[original_file][4],third_it_outfile)
                        protein_file_info[original_file].append(third_it_outfile)
                    else:
                        annotation = sqlite3_search.search_ids_imported(sql_database, "refseq", significant_hits)
                        # Each annotation will have:
                        # query_id gene_id product taxonomy ko_number ko_product
                        # Write the annotations found in the annotation file

                        with open(final_annotation_file, 'a') as final_annotation_fh:
                            for match in annotation:
                                final_annotation_fh.write("{}\t{}\t{}\t{}\t{}\t{}\tNA\tNA\tNA\tNA\tNA\t{}\trefseq\n".format(match[0],
                                match[1], match[2], match[5], match[6], match[3], match[4]))
                                if match[5] != "" and match[5] != "NA":
                                    starting_proteins[str(Path(original_file).name)].remove(match[0])
                                
                        # Check which proteins were not annotated and filter for the next iteration
                        third_it_outfile = str(temporal_protein_folder / (protein_file_info[original_file][0] + ".3it"))
                        if full == True:
                            copyfile(protein_file_info[original_file][4],third_it_outfile)
                            protein_file_info[original_file].append(third_it_outfile)
                        else:
                            # Check if all annotations have ko numbers, if not, get ids to next round of iteration
                            protein_file_info[original_file].append(third_it_outfile)
                            to_retain = starting_proteins[str(Path(original_file).name)].copy()
                            fasta_filter_list.fastA_filter_list(protein_file_info[original_file][2], 
                            third_it_outfile, to_retain, reverse=False)
                            # New structure: 
                            # protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it, Swissprot_Search_File,
                            # filtered_fasta_2it, RefSeq_Search_File, filtered_fasta_3it]
            # When the process is complete, write step completed (+1) in the log folder
            # Also, export dictionary with information to be imported in case of continue  
            process_step += 1
            with open(process_log_folder / "log.txt", 'w') as logfile:
                logfile.write("{}".format(process_step))
            with open(process_log_folder / "structure.pickle", 'wb') as structure_file:
                pickle.dump(protein_file_info, structure_file)
        # --------------------------


        # Finally run remaining proteins against Trembl
    if process_step == 3:
        if light == True:
            print("The light run of MicrobeAnnotator had finished succesfully!")
            print("Summarizing those results...")
            process_step = 5
        else:
            print("Searching proteins against Trembl...\n")
            trembl_database = None
            if method == "blast":
                trembl_database = str(Path(database_folder) / "01.Protein_DB/uniprot_trembl")
            elif method == "diamond":
                trembl_database = str(Path(database_folder) / "01.Protein_DB/uniprot_trembl.dmnd")
            elif method == "sword":
                trembl_database = str(Path(database_folder) / "01.Protein_DB/uniprot_trembl.fasta")
            input_proteins = []
            for infor in protein_file_info.values():
                if len(infor) == 7:
                    input_proteins.append(infor[-1])
            if len(input_proteins) > 0:
                try:
                    pool = multiprocessing.Pool(processes)
                    arguments_to_pass = (output_dir, 'trembl', trembl_database, method,
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
                    if len(info) == 7 and result[0] == info[6]:
                        protein_file_info[filename].append(result[1])
            del temp_dir
            # --------------------------
            # Search annotations in SQLite DB and append to the final annotation file
            print("Extracting Trembl annotation data...\n")
            for original_file, information in protein_file_info.items():
                if len(information) == 8:
                    filtered_sword_results = information[7]
                    final_annotation_file = information[1]
                    significant_hits = []
                    with open(filtered_sword_results) as swissprot_results:
                        for line in swissprot_results:
                            line = line.strip().split()
                            significant_hits.append((line[0], line[1]))
                    if len(significant_hits) == 0:
                        with open(final_annotation_file, 'a') as final_annotation_fh:
                            for original_protein in starting_proteins[str(Path(original_file).name)]:
                                final_annotation_fh.write("{}\tNA\tNo match found\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n".format(original_protein))
                    else:
                        annotation = sqlite3_search.search_ids_imported(sql_database, "trembl", significant_hits)
                        # Each annotation will have:
                        # query_id gene_id product taxonomy ko_number ko_product
                        with open(final_annotation_file, 'a') as final_annotation_fh:
                            for match in annotation:
                                final_annotation_fh.write("{}\t{}\t{}\t{}\tNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\ttrembl\n".format(match[0],
                                match[1], match[3], match[4], match[6], match[7], match[8], match[9], match[10], match[11], match[12]))
                                if match[4] != "":
                                    starting_proteins[str(Path(original_file).name)].remove(match[0])
                        with open(final_annotation_file, 'a') as final_annotation_fh:
                            for original_protein in starting_proteins[str(Path(original_file).name)]:
                                final_annotation_fh.write("{}\tNA\tNo match found\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n".format(original_protein))
            if temporal_protein_folder.is_dir():
                rmtree(temporal_protein_folder)
            # When the process is complete, write step completed (+1) in the log folder
            # Also, export dictionary with information to be imported in case of continue  
            process_step += 1
            with open(process_log_folder / "log.txt", 'w') as logfile:
                logfile.write("{}".format(process_step))
            with open(process_log_folder / "structure.pickle", 'wb') as structure_file:
                pickle.dump(protein_file_info, structure_file)
    # ------------------------

    # Refine annotations with ids from other databases
    if process_step == 4:
        if refine == False:
            process_step += 1
            with open(process_log_folder / "log.txt", 'w') as logfile:
                logfile.write("{}".format(process_step))
            with open(process_log_folder / "structure.pickle", 'wb') as structure_file:
                pickle.dump(protein_file_info, structure_file)
        else:
            import pandas as pd
            print("Improving annotations by searching matches in other databases...")
            # Structure: 
            # protein_file : [protein_file_name, final_annotation_file, filtered_fasta_1it, Swissprot_Search_File,
            # filtered_fasta_2it, RefSeq_Search_File, filtered_fasta_3it, Trembl_Search_File]
            final_annotation_files = []
            for information in protein_file_info.values():
                final_annotation_files.append(str(information[1]))
            
            # Process each annotation file separately
            for annotation_file in final_annotation_files:
                positive_matches = {}
                annotation_table = pd.read_csv(annotation_file, sep="\t", header=0, index_col=0)
                # Parse Uniprot Records
                swissprot_records = annotation_table.loc[(annotation_table['database'] == "swissprot") & (annotation_table['ko_number'] == "NA"), ]
                # Extract records without ko numbers but with ec number
                swissprot_records_ko = swissprot_records.loc[(swissprot_records['ko_number'] == "NA") & (swissprot_records['ec_number'] != "NA"), ]
                swissprot_records_ko_ids = list(swissprot_records_ko['ec_number'])
                if len(swissprot_records_ko_ids) > 0:
                    positive_matches = identifier_conversion.convert_ko_to_ec(swissprot_records_ko_ids, interconversion_database, True)
                    if len(positive_matches) > 0:
                        for identifier, matches in positive_matches.items():
                            annotation_table.loc[(annotation_table['ec_number'] == identifier) & (annotation_table['ko_number'] == "NA"), "ko_number"] = " ".join(matches)
                # Now, extract records without ko number but with interpro ids
                swissprot_records = annotation_table.loc[(annotation_table['database'] == "swissprot") & (annotation_table['ko_number'] == "NA"), ]
                swissprot_records_ko = swissprot_records.loc[(swissprot_records['ko_number'] == "NA") & (swissprot_records['interpro'] != "NA"), ]
                swissprot_records_ko_ids = list(swissprot_records_ko['interpro'])
                if len(swissprot_records_ko_ids) > 0:
                    positive_matches = identifier_conversion.convert_interpro_to_ko(swissprot_records_ko_ids, interconversion_database, False)
                    if len(positive_matches) > 0:
                        for identifier, matches in positive_matches.items():
                            annotation_table.loc[(annotation_table['interpro'] == identifier) & (annotation_table['ko_number'] == "NA"), "ko_number"] = " ".join(matches)
                # Delete used sub-tables
                del swissprot_records
                del swissprot_records_ko
                
                # Parse Trembl Records
                trembl_records = annotation_table.loc[(annotation_table['database'] == "trembl") & (annotation_table['ko_number'] == "NA"), ]
                # Extract records without ko numbers but with ec number
                trembl_records_ko = trembl_records.loc[(trembl_records['ko_number'] == "NA") & (trembl_records['ec_number'] != "NA"), ]
                trembl_records_ko_ids = list(trembl_records_ko['ec_number'])
                if len(trembl_records_ko_ids) > 0:
                    positive_matches = identifier_conversion.convert_ko_to_ec(trembl_records_ko_ids, interconversion_database, True)
                    if len(positive_matches) > 0:
                        for identifier, matches in positive_matches.items():
                            annotation_table.loc[(annotation_table['ec_number'] == identifier) & (annotation_table['ko_number'] == "NA"), "ko_number"] = " ".join(matches)
                # Now, extract records without ko number but with interpro ids
                trembl_records = annotation_table.loc[(annotation_table['database'] == "trembl") & (annotation_table['ko_number'] == "NA"), ]
                trembl_records_ko = trembl_records.loc[(trembl_records['ko_number'] == "NA") & (trembl_records['interpro'] != "NA"), ]
                trembl_records_ko_ids = list(trembl_records_ko['interpro'])
                if len(trembl_records_ko_ids) > 0:
                    positive_matches = identifier_conversion.convert_interpro_to_ko(trembl_records_ko_ids, interconversion_database, False)
                    if len(positive_matches) > 0:
                        for identifier, matches in positive_matches.items():
                            annotation_table.loc[(annotation_table['interpro'] == identifier) & (annotation_table['ko_number'] == "NA"), "ko_number"] = " ".join(matches)
                del trembl_records
                del trembl_records_ko

                # Parse RefSeq Records
                refseq_records = annotation_table.loc[(annotation_table.database == "refseq") & (annotation_table['ko_number'] == "NA"), ]
                # Extract records without ko numbers but with ec number
                refseq_records_ko = refseq_records.loc[(refseq_records['ko_number'] == "NA") & (refseq_records['ec_number'] != "NA"), ]
                refseq_records_ko_ids = list(refseq_records_ko['ec_number'])
                if len(refseq_records_ko_ids) > 0:
                    positive_matches = identifier_conversion.convert_ko_to_ec(refseq_records_ko_ids, interconversion_database, True)
                    if len(positive_matches) > 0:
                        for identifier, matches in positive_matches.items():
                            annotation_table.loc[(annotation_table['ec_number'] == identifier), "ko_number"] = " ".join(matches)
                # Extract refseq ids for records with no ko
                refseq_records = annotation_table.loc[(annotation_table.database == "refseq") & (annotation_table['ko_number'] == "NA"), ]
                refseq_records_ids = list(refseq_records['protein_id'])
                if len(refseq_records_ids) > 0:
                    refseq_records_indices = refseq_records.index
                    query_refseq_pair = list(zip(refseq_records_indices, refseq_records_ids))
                    # Convert refseq to uniprot
                    refseq_uniprot = {}
                    positive_matches = identifier_conversion.convert_refseq_to_uniprot(refseq_records_ids, interconversion_database, False)
                    if len(positive_matches) > 0:
                        for identifier, matches in positive_matches.items():
                            refseq_uniprot[identifier] = matches[0]
                        query_uniprot_pair = []
                        for pair in query_refseq_pair:
                            if pair[1] in refseq_uniprot:
                                query_uniprot_pair.append((pair[refseq_records_indices], refseq_uniprot[pair[1]]))
                        # Now estract the uniprot annotations using the identifiers found in the previous step
                        uniprot_annotations = sqlite3_search.search_ids_imported(sql_database, "trembl", query_uniprot_pair)
                        if len(uniprot_annotations) > 0:
                            # Parse results and modify the original table
                            for annotation in uniprot_annotations:
                                original_index = annotation.pop(0)
                                annotation_table.loc[original_index] = annotation
                # Save the refined table in the original location
                pd.DataFrame.to_csv(annotation_table, sep="\t", index=0, header=True)
            # When the process is complete, write step completed (+1) in the log folder
            # Also, export dictionary with information to be imported in case of continue
            process_step += 1
            with open(process_log_folder / "log.txt", 'w') as logfile:
                logfile.write("{}".format(process_step))
            with open(process_log_folder / "structure.pickle", 'wb') as structure_file:
                pickle.dump(protein_file_info, structure_file)
            print("Done!")
    # ------------------------

    # Parse annotation files and summarize them
    if process_step == 5:
        print("Extracting ko numbers and summarizing results...")
        annotation_files = []
        for information in protein_file_info.values():
            annotation_folder = information[1].parent
            ko_numbers = str(annotation_folder / (information[0] + ".ko"))
            with open(information[1], 'r') as annotations, open(ko_numbers, 'w') as ko_present:
                for line in annotations:
                    line = line.strip().split("\t")
                    # If there are multiple KOs
                    if line[3] != "":
                        if len(line[3].split()) > 1:
                            for element in line[3].split():
                                ko_present.write("{}\n".format(element))
                        else:
                            ko_present.write("{}\n".format(line[3]))
            annotation_files.append(ko_numbers)
        prefix = str(Path(output_dir) / plot_filename)
        regular_modules, bifurcation_modules, structural_modules,  \
        module_information, metabolism_matrix, module_group_matrix = ko_mapper.module_information_importer(annotation_files)
        metabolic_annotation = ko_mapper.global_mapper(regular_modules, bifurcation_modules, structural_modules, annotation_files)
        metabolism_matrix_dropped_relabel, module_colors = ko_mapper.create_output_files(metabolic_annotation, metabolism_matrix, module_information, cluster, prefix)
        ko_mapper.plot_function_barplots(module_colors, module_group_matrix, metabolism_matrix_dropped_relabel, prefix)
        print("MicrobeAnnotator has finished succefully!")
    # ------------------------

if __name__ == "__main__":
    main()

