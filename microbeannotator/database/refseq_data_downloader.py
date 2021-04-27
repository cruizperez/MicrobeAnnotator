#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Date:         Nov 13, 2020

# Description: Download the latest refseq protein dataset from NCBI and
# the accompanying genbank files necessary to create the database.
# ==============================================================================
"""


# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.logging import setup_logger

from shutil import copyfileobj
from shutil import rmtree
from shutil import which
from pathlib import Path
from typing import List
from sys import argv

import multiprocessing
import urllib.request
import subprocess
import argparse
import gzip
import wget
# ==============================================================================


# ==============================================================================
# Initalize logger
# ==============================================================================
logger = setup_logger(__name__)
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
# Function to download refseq data using aspera connect
def refseq_fasta_downloader(
    output_file_folder: Path, ascp_key: str = None) -> None:
    """
    Downloads all RefSeq proteins using Aspera connect
    and merges them in a single output

    Args:
        output_file_folder (Path): Folder to store downloaded protein files.
        ascp_key (str, optional): Aspera SSH key. Defaults to None.
    """
    logger.info("Downloading protein fasta files using Aspera Connect.")
    # Get aspera ssh key
    if ascp_key == None:
        ascp_key = get_aspera_key()
    # Create required folders
    temp_fasta_files = output_file_folder / "temp_refseq_proteins"
    merged_db_folder = output_file_folder / "protein_db"
    temp_fasta_files.mkdir(parents=True, exist_ok=True)
    merged_db_folder.mkdir(parents=True, exist_ok=True)
    # Get latest release for RefSeq and store it
    release_file = urllib.request.urlopen(
        "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER")
    release_number = release_file.readline().decode('utf-8').strip()
    with open(merged_db_folder / "refseq_release.txt", 'w') as relinfo:
        relinfo.write(f"Release number {release_number}")
    # Remove refseq fasta files if present (to avoid repeating files)
    if Path(merged_db_folder / "refseq_protein.fasta").is_file():
        Path(merged_db_folder / "refseq_protein.fasta").unlink()
    # Download compressed protein files
    subprocess.call(["ascp", '-QTr', '-d', '-l', '100M', '-E', '*wgs_mstr*',
                    '-E', '*rna*', '-E', '*genomic*', '-E', '*.gpff.gz',
                    '-i', ascp_key,
                    'anonftp@ftp.ncbi.nlm.nih.gov:/refseq/release/viral/',
                    temp_fasta_files])
    subprocess.call(["ascp", '-QTr', '-d', '-l', '100M', '-E', '*wgs_mstr*',
                    '-E', '*rna*', '-E', '*genomic*', '-E', '*.gpff.gz',
                    '-i', ascp_key,
                    'anonftp@ftp.ncbi.nlm.nih.gov:/refseq/release/bacteria/',
                    temp_fasta_files])
    subprocess.call(["ascp", '-QTr', '-d', '-l', '100M', '-E', '*wgs_mstr*',
                    '-E', '*rna*', '-E', '*genomic*', '-E', '*.gpff.gz',
                    '-i', ascp_key,
                    'anonftp@ftp.ncbi.nlm.nih.gov:/refseq/release/archaea/',
                    temp_fasta_files])
    # Get downloaded files
    viral_list = search_all_files(temp_fasta_files / "viral")
    bacteria_list = search_all_files(temp_fasta_files / "bacteria")
    archaea_list = search_all_files(temp_fasta_files / "archaea")
    final_list = viral_list + bacteria_list + archaea_list
    refseq_proteins = merged_db_folder / "refseq_protein.fasta"
    logger.info("Merging protein files")
    with open(refseq_proteins, 'w') as merged_db:
        for file in final_list:
            with gzip.open(file, 'rt') as temp_file:
                copyfileobj(temp_file, merged_db)
                Path.unlink(file)
    rmtree(temp_fasta_files)
    logger.info("RefSeq proteins successfully downloaded.")

    return refseq_proteins
    
def refseq_fasta_downloader_wget(
    output_file_folder: Path, threads: int) -> None:
    """
    Downloads all RefSeq proteins using wget
    and merges them in a single output

    Args:
        output_file_folder (Path): Folder to store downloaded protein files.
        threads (int): Number of threads to use.
    """
    logger.info("Downloading protein fasta files using wget.")
    # Create required folders
    temp_fasta_files = Path(output_file_folder) / "temp_refseq_proteins"
    merged_db_folder = Path(output_file_folder) / "protein_db"
    Path(temp_fasta_files).mkdir(parents=True, exist_ok=True)
    Path(merged_db_folder).mkdir(parents=True, exist_ok=True)
    # Get latest release for RefSeq and store it
    release_file = urllib.request.urlopen(
        "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER")
    release_number = release_file.readline().decode('utf-8').strip()
    with open(merged_db_folder / "refseq_release.txt", 'w') as relinfo:
        relinfo.write("Release number {}".format(release_number))
    # Remove refseq fasta files if present (to avoid repeating files)
    if Path(merged_db_folder / "refseq_protein.fasta").is_file():
        Path(merged_db_folder / "refseq_protein.fasta").unlink()
    # Download compressed protein files using wget
    file_list = []
    for domain in ["viral", "bacteria", "archaea"]:
        ncbi_url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/" + domain
        domain_url_info = urllib.request.urlopen(ncbi_url)
        for line in domain_url_info:
            line = line.decode('utf-8').strip()
            if "protein.faa.gz" in line:
                filename = str(line.split('"')[1])
                file_url = (
                    f"ftp://ftp.ncbi.nlm.nih.gov/refseq/release/"
                    f"{domain}/{filename}"
                )
                output = f"{str(temp_fasta_files)}/{filename}"
                file_list.append((file_url, output))            
    try:
        pool = multiprocessing.Pool(threads)
        pool.map(refseq_multiprocess_downloader, file_list)
    finally:
        pool.close()
        pool.join()
    # Merge downloaded files
    logger.info("Merging protein files")
    final_list = search_all_files(temp_fasta_files)
    refseq_proteins = merged_db_folder / "refseq_protein.fasta"
    with open(refseq_proteins, 'w') as merged_db:
        for file in final_list:
            with gzip.open(file, 'rt') as temp_file:
                copyfileobj(temp_file,merged_db)
                Path.unlink(file)
    rmtree(temp_fasta_files)
    logger.info("RefSeq proteins successfully downloaded.")

    return refseq_proteins


def refseq_genbank_downloader(
    output_file_folder: Path,
    ascp_key: str = None) -> List[Path]:
    """
    Downloads RefSeq Genbank files with protein metadata using Aspera Connect.

    Args:
        output_file_folder (Path): Folder to store downloaded Genbank files.
        ascp_key (str, optional): Aspera SSH key. Defaults to None.

    Returns:
        List[Path]: List of Genbank files.
    """
    logger.info("Downloading protein genbank files using Aspera Connect.")
    # Get aspera ssh key
    if ascp_key == None:
        ascp_key = get_aspera_key()
    # Make output folder
    temp_gb_dir = output_file_folder / "temp_genbank"
    temp_gb_dir.mkdir(parents=True, exist_ok=True)
    # Download compressed protein files
    subprocess.call(["ascp", '-QTr', '-d', '-l', '100M', '-E', '*wgs_mstr*',
                    '-E', '*rna*', '-E', '*genomic*', '-E', '*.faa.gz',
                    '-i', ascp_key,
                    'anonftp@ftp.ncbi.nlm.nih.gov:/refseq/release/viral/',
                    temp_gb_dir])
    subprocess.call(["ascp", '-QTr', '-d', '-l', '100M', '-E', '*wgs_mstr*',
                    '-E', '*rna*', '-E', '*genomic*', '-E', '*.faa.gz',
                    '-i', ascp_key,
                    'anonftp@ftp.ncbi.nlm.nih.gov:/refseq/release/bacteria/',
                    temp_gb_dir])
    subprocess.call(["ascp", '-QTr', '-d', '-l', '100M', '-E', '*wgs_mstr*',
                    '-E', '*rna*', '-E', '*genomic*', '-E', '*.faa.gz',
                    '-i', ascp_key,
                    'anonftp@ftp.ncbi.nlm.nih.gov:/refseq/release/archaea/',
                    temp_gb_dir])
    # Get downloaded files
    viral_list = search_all_files(temp_gb_dir / "viral")
    bacteria_list = search_all_files(temp_gb_dir / "bacteria")
    archaea_list = search_all_files(temp_gb_dir / "archaea")
    final_list = viral_list + bacteria_list + archaea_list
    logger.info("RefSeq genbank successfully downloaded.")

    return temp_gb_dir, final_list

def refseq_genbank_downloader_wget(
    output_file_folder: Path, threads: int) -> List[Path]:
    """
    Downloads RefSeq Genbank files with protein metadata using wget.

    Args:
        output_file_folder (Path): [description]
        threads (int): [description]

    Returns:
        List[Path]: [description]
    """
    logger.info("Downloading protein genbank files using wget.")
    # Make output folder
    temp_gb_dir = output_file_folder / "temp_genbank"
    temp_gb_dir.mkdir(parents=True, exist_ok=True)
    # Download compressed genbank files using wget
    file_list = []
    for domain in ["viral", "bacteria", "archaea"]:
        ncbi_url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/" + domain
        domain_url_info = urllib.request.urlopen(ncbi_url)
        for line in domain_url_info:
            line = line.decode('utf-8').strip()
            if "protein.gpff.gz" in line:
                filename = str(line.split('"')[1])
                file_url = (
                    f"ftp://ftp.ncbi.nlm.nih.gov/refseq/release/"
                    f"{domain}/{filename}"
                )
                output = f"{str(temp_gb_dir)}/{filename}"
                file_list.append((file_url, output))
    # Download using multiple cores
    try:
        pool = multiprocessing.Pool(threads)
        pool.map(refseq_multiprocess_downloader, file_list)
    finally:
        pool.close()
        pool.join()
    final_list = search_all_files(temp_gb_dir)
    logger.info("RefSeq genbank successfully downloaded.")

    return temp_gb_dir, final_list

def refseq_multiprocess_downloader(file_list: List[str]) -> None:
    """
    Downloads data from refseq using wget

    Args:
        file_list (List[Path]): List of files to download.
    """
    file_url = file_list[0]
    output = file_list[1]
    wget.download(file_url, out=output)

def get_aspera_key() -> str:
    """
    Retrieves Aspera Connect ssh key

    Returns:
        str: Path to Aspera Connect ssh key file
    """
    path_to_ascp = which("ascp")
    install_folder = Path(path_to_ascp).parents[1]
    key = f"{install_folder}/etc/asperaweb_id_dsa.openssh"
    
    return key

def search_all_files(directory: Path) -> List[Path]:
    """
    Search all compressed files in a directory

    Args:
        directory (Path): Directory to search.

    Returns:
        List[Path]: List of files ending in .gz in search directory.
    """
    dirpath = Path(directory)
    assert(dirpath.is_dir())
    file_list = []
    for file in dirpath.iterdir():
        if file.is_file() and file.suffix == ".gz":
            file_list.append(file)
        else:
            continue

    return file_list
# ==============================================================================


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            f"This script downloads the Fasta and Genbank files that\n"
            f"are needed to build the RefSeq annotation database.\n"
            f"By default it downloads both but either can be specified.\n"
            f"Mandatory parameters: -o [output directory]\n"
            f"Optional parameters: See {argv[0]} -h\n"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory")
    mandatory_arguments.add_argument(
        '-o', '--outdir', dest='outdir', action='store', required=True,
        help='Output directory to store the Fasta and Genbank files.')
    # Setup optional arguments
    optional_arguments = parser.add_argument_group("Optional")
    optional_arguments.add_argument(
        '--proteins', dest='proteins', action='store_true', required=False,
        help='Only download the protein Fasta files.')
    optional_arguments.add_argument(
        '--genbank', dest='genbank', action='store_true', required=False,
        help='Only download the Genbank files.')
    optional_arguments.add_argument(
        '--ascp', dest='ascp_key', action='store', required=False, default=None,
        help='Path to Aspera Connect key.')
    # If no arguments are provided
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    # Parse arguments
    outdir = arguments.outdir
    proteins = arguments.proteins
    genbank = arguments.genbank
    ascp_key = arguments.ascp_key

    # Run functions
    if proteins and genbank:
        logger.error(
            f"Please specify only either --proteins or --genbank.\n"
            f"If you want both do not pass any of these flags.")
    if proteins:
        refseq_fasta_downloader(outdir, ascp_key)
    elif genbank:
        refseq_genbank_downloader(outdir, ascp_key)
    else:
        refseq_fasta_downloader(outdir, ascp_key)
        refseq_genbank_downloader(outdir, ascp_key)
# ==============================================================================


# ==============================================================================
# Run main function
# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
