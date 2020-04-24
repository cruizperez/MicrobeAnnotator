#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      0.9
# Date:         March 04, 2020

# Description: Parses compressed dat files and extracts protein
information relevant for annotation purposes.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import gzip

################################################################################
"""---1.0 Define Functions---"""

def parse_uniprot_dat(dat_file, output_table):
    with gzip.open(dat_file, 'rt') as uniprot, open(output_table, 'w') as output:
        gene_id = ""
        accession = ""
        gene_name = ""
        ko_number = ""
        organism = ""
        taxonomy = ""
        function = ""
        compartment = ""
        process = ""
        for line in uniprot:
            if line.startswith("ID", 0):
                gene_id = line.split()[1]
            elif "AC  " in line:
                accession = line.split()[1].replace(";","")
            elif "RecName" in line:
                gene_name = line.split("Full=")[1]
                gene_name = gene_name.split("{")[0].strip().replace(";","")
            elif "OS  " in line:
                organism = ' '.join([organism, line.split("OS")[1].strip()]).replace(".","")
            elif "OC  " in line:
                taxonomy = ' '.join([taxonomy, line.split("OC")[1].strip()]).replace(".","")
            elif "DR   KO;" in line:
                ko_number = line.split()[2].replace(";", "")
            elif "DR   GO;" in line:
                if "; F:" in line:
                    function = ''.join([function, line.split("GO;")[1].strip(), " -- "])
                elif "; C:" in line:
                    compartment = ''.join([compartment, line.split("GO;")[1].strip(), " -- "])
                elif "; P:" in line:
                    process = ''.join([process, line.split("GO;")[1].strip(), " -- "])
            elif "//\n" in line:
                output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene_id, 
                accession, gene_name, ko_number, organism, taxonomy, function, compartment, process))
                gene_id = ""
                accession = ""
                gene_name = ""
                ko_number = ""
                organism = ""
                taxonomy = ""
                function = ""
                compartment = ""
                process = ""

################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''This script parses a Uniprot.dat file and outputs a table with\n'''
                                                    '''the ID, Accession, Gene Name, Organism, Taxonomy, KEGG ID, Function, Compartment, and Process.\n
                                                    For faster usage in alrge files use gnu parallel (read script file to see how)\n'''
                                    '''\nGlobal mandatory parameters: [Input Uniprot.dat File]\n'''
                                    '''\nOptional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_dat', action='store', required=True, help='Uniprot.dat file to parse')
    parser.add_argument('-o', '--output', dest='output_table', action='store', required=False, help='Output table')
    args = parser.parse_args()

    input_dat = args.input_dat
    output_table = args.output_table

    parse_uniprot_dat(input_dat, output_table)

if __name__ == "__main__":
    main()
