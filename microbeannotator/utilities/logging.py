#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz-Perez
# Email:        cruizperez3@gatech.edu
# Institution:  Georgia Institute of Technology
# Date:         11 April 2021

# Description: This script defines the logger and its options
# ==============================================================================
"""

# ==============================================================================
# Import modules
# ==============================================================================
import logging
import socket
# ==============================================================================


# ==============================================================================
# Define logging methods and filter classes
# ==============================================================================
# Class to add filters for handlers
class LoggerFilter(object):
    def __init__(self, level):
        self.__level = level
    
    def filter(self, logRecord):
        return logRecord.levelno <= self.__level

# Setup logger with handlers for different levels
def setup_logger(name, level: str = 'INFO'):
    # Define format for handlers
    info_format = (f"%(asctime)s " + "[%(levelname)s]: %(message)s")
    warning_format = (
        f"%(asctime)s {socket.gethostname()} %(name)s: "
        f"Function: %(funcName)s - Line: %(lineno)s "
        + "[%(levelname)s]: %(message)s")
    # Define handler for information
    info_handler = logging.StreamHandler()
    info_handler.setLevel(logging.INFO)
    info_handler.setFormatter(logging.Formatter(info_format))
    info_handler.addFilter(LoggerFilter(logging.INFO))
    # Define handler for warnings
    warning_handler = logging.StreamHandler()
    warning_handler.setLevel(logging.WARNING)
    warning_handler.setFormatter(logging.Formatter(warning_format))
    warning_handler.addFilter(LoggerFilter(logging.WARNING))
    # Define handler for errors and critical
    error_handler = logging.StreamHandler()
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(logging.Formatter(warning_format))
    # Initialize logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(info_handler)
    logger.addHandler(warning_handler)
    logger.addHandler(error_handler)

    return logger
# ==============================================================================