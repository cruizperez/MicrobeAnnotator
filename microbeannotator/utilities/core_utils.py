#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz-Perez
# Email:        cruizperez3@gatech.edu
# Institution:  Georgia Institute of Technology
# Date:         11 April 2021

# Description: This script contains utilities used by scripts in several
# locations.
# ==============================================================================
"""

# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.logging import setup_logger

from random import choices
from shutil import which
from pathlib import Path
from typing import List
from sys import exit

import psutil
import string
# ==============================================================================


# ==============================================================================
# Initalize logger
# ==============================================================================
logger = setup_logger(__name__)
# ==============================================================================


# ==============================================================================
# Define functions 
# ==============================================================================
# Generate random strings
def random_string(length: int) -> str:
    """
    Generate random string

    Args:
        length (int): Lenght of the string.

    Returns:
        str: Random string of length [lenght]
    """
    return ''.join(choices(string.ascii_uppercase + string.digits, k = length))

# Validate input parameters
def input_validator(
    method: str, method_bin: Path, input_list: List[str], file_list: Path,
    cluster: str, processes: int, threads: int):
    # Validate method
    method_options = ['blast', 'sword', 'diamond']
    if method not in method_options:
        logger.error(
            f"Invalid search method: {method}.\n"
            f"Please select one of 'diamond' 'blast', or 'sword'")
        exit(1)
    # Validate method call
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
        method_call = str(method_bin / call_name)
    if which(method_call) == None:
            if method_bin == None:
                logger.error(
                    f"{method_call} not found in PATH.\n"
                    f"Please provide the correct path to the "
                    f"directory with {method_call} binaries.")
                exit(1)
            else:
                logger.error(
                    f"{method_call} not found in {method_bin}.\n"
                    f"Please provide the correct path to the "
                    f"directory with {method_call} binaries.")
                exit(1)
    # Validate input files
    if input_list != None and file_list != None:
        logger.error(
            f"Please provide only a list of files (--input) or an\n"
            f"input list (--list), not both.\n")
        exit(1)
    elif input_list != None:
        for file in input_list:
            if not Path(file).is_file():
                logger.error(
                    f"File {file} does not exist. Please check"
                    f"the file exists and can be accessed."
                )
                exit(1)
    elif file_list != None:
        with open(file_list, 'r') as infile:
            for line in infile:
                if not Path(line.strip()).is_file():
                    logger.error(
                        f"File {line.strip()} does not exist. Please check"
                        f"the file exists and can be accessed."
                    )
                    exit(1)
    # Validate number of processes and threads
    max_cores = psutil.cpu_count()
    if (processes * threads) > max_cores:
        logger.warning(
            f"The number of processes {processes} and threads\n"
            f"{threads} exceeds the maximum number of cores\n"
            f"available in the system {max_cores}. This can result\n"
            f"in a decrease in performance.")
    # Validate cluster request
    cluster_options = ['cols', 'rows', 'both']
    if cluster != None and cluster not in cluster_options:
        logger.error(
            f"Invalid cluster option: {cluster}.\n"
            f"Please select one of 'cols' 'rows', or 'both'\n"
            f"or don't pass this flag for no clustering.")
        exit(1)
    
# ==============================================================================