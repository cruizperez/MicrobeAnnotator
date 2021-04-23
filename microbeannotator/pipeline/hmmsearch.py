#!/usr/bin/env python

"""
# ==============================================================================
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.9
# Date:		   14 February 2020

# Description: This script runs hmmsearch and filters the results based on the
# thresholds described by kofamscan.
# ==============================================================================
"""

# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.core_utils import random_string
from microbeannotator.utilities.logging import setup_logger
from microbeannotator.database import kofam_profile_parser as kofam

from tempfile import NamedTemporaryFile
from tempfile import TemporaryDirectory
from pathlib import Path
from typing import List
from typing import Tuple

import multiprocessing
import itertools
import hmmer
# ==============================================================================


# ==============================================================================
# Initalize logger
# ==============================================================================
logger = setup_logger(__name__)
# ==============================================================================


# ==============================================================================
# Define functions 
# ==============================================================================
def hmmer_search(protein_model: Tuple[Path, Path]) -> List[hmmer.tbl.TBLRow]:
    # Get protein and model paths
    protein_path = protein_model[0]
    model_path = protein_model[1]
    # Create model database
    model_db = hmmer.HMMER(model_path)
    # Create temporary directory
    temp_dir = TemporaryDirectory(suffix='_microbeannotator')
    temp_file = NamedTemporaryFile(suffix=random_string(20))
    temporal_result = Path(temp_dir.name) / Path(temp_file.name)
    search_result = model_db.search(
        protein_path,
        output=temporal_result)
    if len(search_result.tbl) > 0:
        search_result_filt = hmmer_filter(search_result.tbl)
    
        return search_result_filt


def hmmer_filter(
    hmmsearch_result: List[hmmer.tbl.TBLRow]) -> List[hmmer.tbl.TBLRow]:
    # Initialize list with filtered results
    hmmsearch_result_filt = []
    # Get model name, threshold and score_type
    model_name = hmmsearch_result[0].query.name
    threshold = hmm_model_info[model_name].threshold
    score_type = hmm_model_info[model_name].score_type
    if score_type == 'full':
        for result in hmmsearch_result:
            if float(result.full_sequence.score) >= threshold:
                hmmsearch_result_filt.append(result)
    elif score_type == 'domain':
        for result in hmmsearch_result:
            if float(result.best_1_domain.score) >= threshold:
                hmmsearch_result_filt.append(result)
    else:
        for result in hmmsearch_result:
            hmmsearch_result_filt.append(result)

    return hmmsearch_result_filt


def write_hmmsearch_results(
    combined_hmmsearch_results: List[hmmer.tbl.TBLRow],
    out_path: Path) -> None:
    with open(out_path, 'w') as output:
        output.write((
            "#\t\t\t\t\t----\tfull_sequence\t----\t----\tbest_1_domain\t----\t"
            "---\t---\t---\tdomain\tnumber\t"
            "estimation\t---\t---\tdescription_of_target\n"
            "# query_name\taccession\ttarget_name\taccession\tscore_type\t"
            "E-value\tscore\tbias\tE-value\tscore\tbias\texp\treg\tclu\tov\t"
            "env\tdom\trep\tinc\tdescription_of_target\n"
            "# -----------\t---------\t----------\t---------\t----------\t"
            "-------\t-----\t----\t-------\t-----\t----\t---\t---\t---\t--"
            "\t---\t---\t---\t---\t---------------------\n"))
        for result in combined_hmmsearch_results:
            if result == None:
                continue
            elif len(result) == 0:
                continue
            else:
                for record in result:
                    output.write((
                        f"{record.target.name}\t{record.target.accession}\t"
                        f"{record.query.name}\t{record.query.accession}\t"
                        f"{hmm_model_info[record.query.name].score_type}\t"
                        f"{record.full_sequence.e_value}\t"
                        f"{record.full_sequence.score}\t"
                        f"{record.full_sequence.bias}\t"
                        f"{record.best_1_domain.e_value}\t"
                        f"{record.best_1_domain.score}\t"
                        f"{record.best_1_domain.bias}\t"
                        f"{record.domain_numbers.exp}\t"
                        f"{record.domain_numbers.reg}\t"
                        f"{record.domain_numbers.clu}\t"
                        f"{record.domain_numbers.ov}\t"
                        f"{record.domain_numbers.env}\t"
                        f"{record.domain_numbers.dom}\t"
                        f"{record.domain_numbers.rep}\t"
                        f"{record.domain_numbers.inc}\t"
                        f"{hmm_model_info[record.query.name].definition}\n"))


def best_match_selector(raw_results: Path) -> Path:
    outfile = raw_results.parent / (raw_results.name + '.filt')
    with open(raw_results, 'r') as infile, \
        open(outfile, 'w') as output:
        best_matches = {}
        for line in infile:
            if line.startswith('#'):
                output.write(line)
            else:
                record = line.strip().split('\t')
                if record[0] not in best_matches:
                    best_matches[record[0]] = [record[6], line]
                else:
                    if record[6] > best_matches[record[0]][0]:
                        best_matches[record[0]] = [record[6], line]
        for record in best_matches.values():
            output.write(record[1])
    return outfile


def protein_to_model_mapper(
    protein_list: List[Path],
    model_list: List[Path]) -> List[Tuple[Path, Path]]:
    combined_list = []
    for path_product in itertools.product(*[protein_list, model_list]):
        combined_list.append(path_product)
    
    return combined_list


def write_hmmsearch_annotation(
    hmmsearch_results: Path, annotation_file: Path):
    with open(hmmsearch_results) as infile, \
        open(annotation_file, 'a') as output:
        for line in infile:
            if line.startswith("#"):
                continue
            else:
                line = line.strip().split('\t')
                query_id = line[0]
                ko_id = line[2]
                ko_product = line[19]
                output.write(
                    f"{query_id}\tNA\tNA\t{ko_id}\t{ko_product}\tNA\t"
                    f"NA\tNA\tNA\tNA\tNA\tNA\tkofam\n"
                )


def hmmsearch_launcher(
    protein_model_list: List[Tuple[Path, Path]],
    ko_list_db: Path,
    output_path: Path,
    processes: int,
    threads: int):
    # Launch hmmsearch in one process per protein file vs model
    global hmm_model_info
    try:
        hmm_model_info
    except NameError:
        hmm_model_info = kofam.parse_profiles(ko_list_db)
    try:
        pool_search = multiprocessing.Pool(processes * threads)
        combined_hmmsearch_results = pool_search.map(
            hmmer_search, protein_model_list)
    finally:
        pool_search.close()
        pool_search.join()
    
    logger.info("Filtering and saving results...")
    # Get results and write them to output file
    write_hmmsearch_results(combined_hmmsearch_results, output_path)
    # Read intermediate results and select best match per protein
    best_match_file = best_match_selector(output_path)
    # Get annotated and hypothetical proteins
    annotated_proteins = []
    hypothetical_proteins = []
    with open(best_match_file) as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            else:
                line = line.strip().split('\t')
                if "hypothetical" not in line[19].lower() and \
                    "uncharacterized" not in line[19].lower() and \
                    "unknown function" not in line[10].lower():
                    annotated_proteins.append(line[0])
                else:
                    hypothetical_proteins.append(line[0])
    logger.info("Search against KOfam done!")
    return best_match_file, annotated_proteins, hypothetical_proteins
