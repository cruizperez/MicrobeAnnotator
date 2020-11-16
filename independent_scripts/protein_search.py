#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Iteratively annotates proteins using different databases
# and multiple search methods.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import wget
from shutil import copyfileobj
from pathlib import Path
import urllib.request
import gzip
from shutil import which
from shutil import rmtree

################################################################################
"""---1.0 Define Functions---"""
def kofamscan_annotation(protein_file, multi_argument):
    # Import modules
    import subprocess
    import random
    import string
    # Unpack arguments
    output_directory = multi_argument[0]
    threads = multi_argument[1]
    bin_path = multi_argument[2]
    # Create folders
    output_directory = Path(output_directory)
    kofamscan_result_folder = output_directory / 'kofamscan_results'
    kofamscan_result_folder.mkdir(parents=True, exist_ok=True)
    protein_file_name = Path(protein_file).name
    # Call kofamscan and filter results
    if bin_path == None:
        kofam_call = 'exec_annotation'
    else:
        kofam_call = str(Path(bin_path) / 'exec_annotation')
    kofam_output_file = str(kofamscan_result_folder / (protein_file_name + '.kofam'))
    temp_folder = ''.join([random.choice(string.ascii_letters + string.digits) for n in range(15)])
    subprocess.call([kofam_call, '-o', kofam_output_file, '--cpu', str(threads),
                    '--tmp-dir', temp_folder, protein_file])
    filtered_output_file = kofam_output_file + '.filt'
    kofamscan_filter(kofam_output_file, filtered_output_file)
    if Path(temp_folder).is_dir():
        rmtree(temp_folder)

    # Extract ids of proteins annotated and add annotations into final file
    final_annotation_folder = output_directory / 'annotation_results'
    final_annotation_folder.mkdir(parents=True, exist_ok=True)
    final_file = Path(final_annotation_folder) / (protein_file_name + '.annotations')
    ids_proteins_annotated = []
    #* The final annotation will have the following columns
    #* query_id protein_id product ko_number ko_product taxonomy function compartment process EC InterPro Pfam database
    with open(filtered_output_file, 'r') as ko_annotated, open(final_file, 'w') as output:
        output.write("query_id\tprotein_id\tproduct\tko_number\tko_product\ttaxonomy\tfunction_go\tcompartment_go\tprocess_go\tinterpro\tpfam\tec_number\tdatabase\n")
        for line in ko_annotated:
            if "#" in line:
                continue
            else:
                line = line.strip().split("\t")
                ids_proteins_annotated.append(line[1])
                output.write("{}\tNA\tNA\t{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tkofamscan\n".format(line[1], line[2], line[6]))

    return (protein_file, protein_file_name, ids_proteins_annotated, final_file)


def kofamscan_filter(kofamscan_input, outfile):
    # Import modules
    from random import randint
    from sys import exit

    header_list = []
    kofamscan_results = {}
    with open(kofamscan_input, 'r') as input:
        for line in input:
            if line.startswith("#"):
                if 'gene name' in line:
                    headers = []
                    headers.append('')
                    head = line.strip().split()
                    name = ' '.join(head[0:3])
                    headers.append(name)
                    headers += head[3:]
                    header_list.append("\t".join(headers))
                else:
                    header = ['']
                    header += line.strip().split()
                    header_list.append("\t".join(header))
            elif line.startswith('*'):
                result = line.strip().split()
                gene_name = result[1]
                annot_desc = ' '.join(result[6:])
                if gene_name not in kofamscan_results:
                    kofamscan_results[gene_name] = result[0:6]
                    kofamscan_results[gene_name].append(annot_desc)
                else:
                    if result[4] > kofamscan_results[gene_name][4]:
                        kofamscan_results[gene_name] = result[0:6]
                        kofamscan_results[gene_name].append(annot_desc)
                    elif result[4] == kofamscan_results[gene_name][4]:
                        if randint(0,1) > 0:
                            kofamscan_results[gene_name] = result[0:6]
                            kofamscan_results[gene_name].append(annot_desc)
                        else:
                            continue
                    else:
                        continue
            else:
                continue

    with open(outfile, 'w') as output_file:
        for element in header_list:
            output_file.write("{}\n".format(element))
        for hit in kofamscan_results.values():
            output_file.write("{}\n".format("\t".join(hit)))


def blast_filter_slow(input_tab, outfile, id_perc, bitscore, evalue, aln_percent):
    # Import modules
    from random import choice

    blast_hits = {}
    # Retrieve best matches
    with open(input_tab) as blast_input:
        for line in blast_input:
            line = line.strip()
            hit = line.split("\t")
            if float(hit[2]) >= id_perc and float(hit[11]) >= bitscore \
            and float(hit[10]) <= evalue and float(hit[3])*100/float(hit[12]) >= aln_percent:
                if hit[0] not in blast_hits:
                    blast_hits[hit[0]] = [float(hit[11]), [line]]
                else:
                    if float(hit[11]) < blast_hits[hit[0]][0]:
                        continue
                    elif float(hit[11]) > blast_hits[hit[0]][0]:
                        blast_hits[hit[0]] = [float(hit[11]), [line]]
                    else:
                        blast_hits[hit[0]][1].append(line)
            else:
                continue
    with open(outfile, 'w') as output:
        for hit_values in blast_hits.values():
            output.write("{}\n".format(choice(hit_values[1])))


def add_query_len(input_tab, protein_file, outfile):
    # Import modules
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    # Read protein lenghts
    protein_len = {}
    with open(protein_file, 'r') as proteins:
        for title, seq in SimpleFastaParser(proteins):
            protein_len[title.split()[0]] = len(seq)
    with open(input_tab, 'r') as blast_like, open(outfile, 'w') as output:
        for line in blast_like:
            line = line.strip()
            protein = line.split("\t")[0]
            output.write("{}\t{}\n".format(line, protein_len[protein]))

def similarity_search(protein_file, multiple_arguments):
    # Import modules
    import subprocess
    from shutil import which
    from sys import exit
    # Unpack arguments
    output_directory = multiple_arguments[0]
    db_version = multiple_arguments[1]
    database = multiple_arguments[2]
    method = multiple_arguments[3]
    threads = multiple_arguments[4]
    id_perc = multiple_arguments[5]
    bitscore = multiple_arguments[6]
    evalue = multiple_arguments[7]
    aln_percent = multiple_arguments[8]
    bin_path = multiple_arguments[9]
    # Create folders
    output_directory = Path(output_directory)
    db_result_folder = output_directory / (db_version + '_results')
    db_result_folder.mkdir(parents=True, exist_ok=True)
    protein_file_name = Path(protein_file).name

    # Call method and filter results
    if method == "blast":
        if bin_path == None:
            blast_call = "blastp"
        else:
            blast_call = str(Path(bin_path) / "blastp")
        if which(blast_call) == None:
            if bin_path == None:
                exit("Blastp not found in PATH, please provide the correct path to the folder with binaries")
            else:
                exit("Blastp not found in " + bin_path + ", please provide the correct path to the folder with binaries")
        else:
            method_output_file = str(db_result_folder / (protein_file_name + '.' + method))
            subprocess.call([blast_call, "-query", protein_file, "-db", database, "-out", method_output_file,
            "-max_target_seqs", "6", "-num_threads", str(threads), "-outfmt", """6 std qlen slen"""])
            filtered_file = method_output_file + '.filt'
            blast_filter_slow(method_output_file, filtered_file, id_perc, bitscore, evalue, aln_percent)
            return (protein_file, filtered_file)
    elif method == 'diamond':
        if bin_path == None:
            diamond_call = "diamond"
        else:
            diamond_call = str(Path(bin_path) / "diamond")
        # if which(method) == None:
        #     if bin_path == None:
        #         exit("Diamond not found in PATH, please provide the correct path to the folder with binaries")
        #     else:
        #         exit("Diamond not found in " + bin_path + ", please provide the correct path to the folder with binaries")
        # else:
        method_output_file = str(db_result_folder / (protein_file_name + '.' + method))
        subprocess.call([diamond_call, 'blastp', '--db', database, '--out', method_output_file, 
        '--outfmt', "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore", "qlen", "slen", "--threads", str(threads),
        "--unal", "0", "--max-target-seqs", "6", "--query", protein_file])
        filtered_file = method_output_file + '.filt'
        blast_filter_slow(method_output_file, filtered_file, id_perc, bitscore, evalue, aln_percent)
        return (protein_file, filtered_file)
    elif method == 'sword':
        if bin_path == None:
            sword_call = "sword"
        else:
            sword_call = str(Path(bin_path) / "sword")
        if which(sword_call) == None:
            if bin_path == None:
                exit("Sword not found in PATH, please provide the correct path to the folder with binaries")
            else:
                exit("Sword not found in " + bin_path + ", please provide the correct path to the folder with binaries")
        else:
            method_output_file = str(db_result_folder / (protein_file_name + '.' + method))
            subprocess.call([sword_call, '-i', protein_file, '-j', database, 
            '-o', method_output_file, "-f", "bm8", "-a", "10000", "-k", "5", "-c", "15000",
            "-T", "16", "-t", str(threads)])
            temporal_file = method_output_file + '.temp'
            add_query_len(method_output_file, protein_file, temporal_file)
            filtered_file = method_output_file + '.filt'
            blast_filter_slow(temporal_file, filtered_file, id_perc, bitscore, evalue, aln_percent)
            Path(temporal_file).unlink()
            return (protein_file, filtered_file)






