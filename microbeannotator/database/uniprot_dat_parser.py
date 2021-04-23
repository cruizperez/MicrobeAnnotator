#!/usr/bin/env python

"""
# ==============================================================================
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Parses compressed dat files and extracts protein
information relevant for annotation purposes.
# ==============================================================================
"""

# ==============================================================================
# Import modules
# ==============================================================================
from microbeannotator.utilities.logging import setup_logger

from pathlib import Path
from sys import argv

import argparse
import gzip
# ==============================================================================


# ==============================================================================
# Initalize logger
# ==============================================================================
logger = setup_logger(__name__)
# ==============================================================================


# ==============================================================================
# Define functions
# ==============================================================================
# Function to parse UniProt dat files
def parse_uniprot_dat(dat_file: Path, output_file_table: Path):
    output_folder = output_file_table.parent
    uniprot_to_refseq = Path(output_folder) / "uniprot_to_refseq.txt"
    logger.info(f"Parsing {dat_file} and storing in {output_file_table}")
    with gzip.open(dat_file, 'rt') as uniprot, \
        open(output_file_table, 'w') as output_file, \
        open(uniprot_to_refseq, 'a') as uni_to_ref:
        gene_id = ""
        accession = ""
        gene_name = ""
        ko_number = ""
        organism = ""
        taxonomy = ""
        function = ""
        compartment = ""
        process = ""
        interpro = ""
        pfam = ""
        ec_number = ""
        refseq_code = ""
        for line in uniprot:
            if line.startswith("ID", 0):
                gene_id = line.split()[1]
            elif "AC  " in line:
                accession = line.split()[1].replace(";","")
            elif line.startswith("DE") and "Full=" in line:
                gene_name = line.split("Full=")[1]
                gene_name = gene_name.split("{")[0].strip().replace(";","")
                gene_name = gene_name.lower()
            elif "OS  " in line:
                line = line.split("OS")[1].strip()
                organism = ' '.join([organism, line]).replace(".","")
            elif "OC  " in line:
                line = line.split("OC")[1].strip()
                taxonomy = ' '.join([taxonomy, line]).replace(".","")
            elif "DR   KO;" in line:
                ko_number = line.split()[2].replace(";", "")
            elif "DR   GO;" in line:
                if "; F:" in line:
                    code = line.strip().split(";")[1]
                    code = code.strip()
                    if function == "":
                        function = code
                    else:
                        function = ''.join([function, " ", code])
                elif "; C:" in line:
                    code = line.strip().split(";")[1]
                    code = code.strip()
                    if compartment == "":
                        compartment = code
                    else:
                        compartment = ''.join([compartment, " ", code])
                elif "; P:" in line:
                    code = line.strip().split(";")[1]
                    code = code.strip()
                    if process == "":
                        process = code
                    else:
                        process = ''.join([process, " ", code])
            elif "DR   InterPro" in line:
                code = line.strip().split()[2]
                code = code.replace(";", "")
                if interpro == "":
                    interpro = code
                else:
                    interpro = ' '.join([interpro, code])
            elif "DR   Pfam" in line:
                code = line.strip().split()[2]
                code = code.replace(";", "")
                if pfam == "":
                    pfam = code
                else:
                    pfam = ' '.join([pfam, code])
            elif line.startswith("DE") and "EC=" in line:
                ec_code = line.strip().split()[1]
                ec_code = ec_code.replace(";", "")
                ec_code = ec_code.replace("EC=", "")
                if ec_number == "":
                    ec_number = ec_code
                else:
                    ec_number = ' '.join([ec_number, ec_code])
            elif line.startswith("DR") and "RefSeq" in line:
                refseq_code = line.strip().split()[2]
                refseq_code = refseq_code.replace(";", "")
                uni_to_ref.write("{}\t{}\n".format(gene_id, refseq_code))
            elif "//\n" in line:
                if ko_number == "":
                    ko_number = "NA"
                if organism == "":
                    organism = "NA"
                if function == "":
                    function = "NA"
                if compartment == "":
                    compartment = "NA"
                if process == "":
                    process = "NA"
                if interpro == "":
                    interpro = "NA"
                if pfam == "":
                    pfam = "NA"
                if ec_number == "NA":
                    ec_number = "NA"
                output_file.write(
                    f"{gene_id}\t{accession}\t{gene_name}\t{ko_number}\t"
                    f"{organism}\t{taxonomy}\t{function}\t{compartment}\t"
                    f"{process}\t{interpro}\t{pfam}\t{ec_number}\n")
                gene_id = ""
                accession = ""
                gene_name = ""
                ko_number = ""
                organism = ""
                taxonomy = ""
                function = ""
                compartment = ""
                process = ""
                interpro = ""
                pfam =""
                ec_number = ""
                refseq_code = ""
    logger.info(f"Finished")

    return uniprot_to_refseq
# ==============================================================================


# ==============================================================================
# Define main function
# ==============================================================================
def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            f"This script parses a UniProt.dat file and returns a \n"
            f"tab-separated table with: \n"
            f"GeneID - Accession - Gene Name - KO Number - Organism - "
            f"Taxonomy - Function GO - Compartment GO - Process GO - "
            f"InterPro ID - Pfam ID - EC Number\n"
            f"Mandatory parameters: -i [input .dat file] -o [output table]\n"
            f"Optional parameters: See {argv[0]} -h\n"))
    # Setup mandatory arguments
    mandatory_arguments = parser.add_argument_group("Mandatory")
    mandatory_arguments.add_argument(
        '-i', '--input', dest='input_dat', action='store',
        required=True, help='Uniprot.dat file to parse')
    mandatory_arguments.add_argument(
        '-o', '--output_file', dest='output_file_table',action='store',
        required=True, help='output_file table')
    # If no arguments are provided
    if len(argv) == 1:
        parser.print_help()
        exit(0)
    arguments = parser.parse_args()

    # Parse arguments
    input_dat = arguments.input_dat
    output_file_table = arguments.output_file_table
    # Run functions
    parse_uniprot_dat(input_dat, output_file_table)
# ==============================================================================


# ==============================================================================
# Run main function
# ==============================================================================
if __name__ == "__main__":
    main()
# ==============================================================================
