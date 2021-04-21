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
import string
from random import choices
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
# ==============================================================================