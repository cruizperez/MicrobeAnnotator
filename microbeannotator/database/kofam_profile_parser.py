#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz-Perez
# Email:        cruizperez3@gatech.edu
# Institution:  Georgia Institute of Technology
# Date:         11 April 2021

# Description: This script downloads the profiles and scores for all hmm models
# used in KOfamscan
# ==============================================================================
"""

# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.logging import setup_logger

from collections import namedtuple
from collections import Counter
from datetime import date
from typing import Tuple
from typing import Dict
from pathlib import Path
from sys import argv

import subprocess
import argparse
import shutil
import wget
# ==============================================================================


# ==============================================================================
# TODO: Upload data to our server for backup
# ==============================================================================


# ==============================================================================
# Initalize logger
# ==============================================================================
logger = setup_logger(__name__)
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
# Download latest KOfamscan profiles and scores
def kofam_downloader(outdir: Path) -> Tuple[Path, Path]:
    """
    Downloads and extracts KOfam data

    Args:
        outdir (Path): Directory of database to store data.

    Returns:
        Tuple[Path, Path]: Location of processed data
    """
    logger.info("Downloading KOfam data")
    # Set download paths
    db_dir = outdir / 'kofam_data'
    db_dir.mkdir(parents=True, exist_ok=True)
    profile_info_path = db_dir / 'ko_list.gz'
    ko_profiles_path = db_dir / 'profiles.tar.gz'
    profile_info_path_extr = db_dir / 'ko_list'
    ko_profiles_path_extr = db_dir / 'profiles'
    # Download latest annotations
    try:
        wget.download("ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz",
                            out=str(profile_info_path))
        subprocess.call(['gunzip', str(profile_info_path)])
    except Exception:
        logger.warning(
            f"Data not available from kofam source. "
            f"Using backup data from XXX")
    # Download latest profiles
    try:
        wget.download("ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz",
                            out=str(ko_profiles_path))
        logger.info("Extracting profiles")
        subprocess.call(['tar', 'xfz', str(ko_profiles_path),
        '-C', str(db_dir)])
        ko_profiles_path.unlink()
    except Exception:
        logger.warning(
            f"Data not available from kofam source. "
            f"Using backup data from XXX")
    
    # Write date accessed
    today = date.today()
    download_date = outdir / 'kofam_data/date_accessed.txt'
    with open(download_date, 'w') as output:
        output.write(f'{today.strftime("%B %d, %Y")}')
    logger.info("Finished")

    return profile_info_path_extr, ko_profiles_path_extr

# Format hmm models into groups
def kofam_formatter(ko_profiles_path_extr: Path) -> None:
    """
    Formats data downloaded for kofamscan

    Args:
        ko_profiles_path_extr (Path): Path to location of hmm profiles.
    """
    logger.info("Formatting KOfam profiles")  
    prokaryote_list = []
    eukaryote_list = []
    common_list = []
    with open(ko_profiles_path_extr / 'prokaryote.hal', 'r') as infile:
        for line in infile:
            prokaryote_list.append(line.strip())
    with open(ko_profiles_path_extr / 'eukaryote.hal', 'r') as infile:
        for line in infile:
            eukaryote_list.append(line.strip())
    # Get models shared by prokaryote and eukaryote
    for model, occurence in dict(Counter(prokaryote_list + eukaryote_list)).items():
        if occurence > 1:
            common_list.append(model)
    # Get unique modules per prokaryote or eukaryote
    prokaryote_only = list(set(prokaryote_list).difference(set(common_list)))
    eukaryote_only = list(set(eukaryote_list).difference(set(common_list)))
    # Get models not in a subset
    accounted_models = eukaryote_only + prokaryote_only + common_list
    total_models = []
    for file in ko_profiles_path_extr.iterdir():
        if file.suffix == '.hmm':
            total_models.append(file.name)
    orfan_models = list(set(total_models).difference(set(accounted_models)))
    # Write common models into groups of 500
    hmm_datasets = [common_list, prokaryote_only, eukaryote_only, orfan_models]
    hmm_datasets_names = ['common', 'prokaryote', 'eukaryote', 'independent']
    chunk_size = 500
    for hmm_list, name in zip(hmm_datasets, hmm_datasets_names):
        file_counter = 1
        file_list = []
        for i in range(0, len(hmm_list), chunk_size):
            group = hmm_list[i:i+chunk_size]
            with open(ko_profiles_path_extr/ f"{name}_{file_counter}.model", 'w') \
                as output:
                for file in group:
                    with open(ko_profiles_path_extr / file, 'r') as temporal:
                        shutil.copyfileobj(temporal, output)
                    (ko_profiles_path_extr / file).unlink()
            file_list.append(f"{name}_{file_counter}.model")
            file_counter += 1
        with open(ko_profiles_path_extr / f"{name}.list", 'w') as output:
            for file in file_list:
                output.write(f"{str(file)}\n")
    logger.info("Finished") 


# Parse kofamscan model info
def parse_profiles(profile_info_path: Path) -> Dict[str, namedtuple]:
    """
    Parses profile thresholds, score_type and definitions from kofamscan
    ko_list.

    Args:
        profile_info (Path): Path to ko_list file with profile information.

    Returns:
        Dict[str, namedtuple]: Dictionary with profile information.
    """
    logger.info("Parsing HMM profile metadata") 
    profile_data: Dict[str, namedtuple] = {}
    Profile = namedtuple('Profile', ["threshold", "score_type", "definition"])
    with open(profile_info_path, 'r') as infile:
        for line in infile:
            if line.startswith("knum"):
                continue
            else:
                line = line.strip().split("\t")
                if line[1] == '-':
                    profile_data[line[0]] = Profile(
                        line[1],line[2], line[11])
                else:
                    profile_data[line[0]] = Profile(
                        float(line[1]),line[2], line[11])

    return profile_data
# ==============================================================================


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            f"This script downloads and parses KOfam HMM profiles\n"
            f"to be used by MicrobeAnnotator.\n"
            f"Mandatory parameters: -o [output directory]\n"))
    # Setup argument group
    mandatory_arguments = parser.add_argument_group("Mandatory")
    mandatory_arguments.add_argument(
        '-o', '--outdir', dest='outdir',action='store',
        required=True, help='Output directory.')
    # If no arguments are provided
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    # Parse arguments
    outdir = arguments.outdir
    # Run functions
    profile_info_path_extr, ko_profiles_path_extr = kofam_downloader(outdir)
    kofam_formatter(ko_profiles_path_extr)
# ==============================================================================


# ==============================================================================
# Run main function
# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================