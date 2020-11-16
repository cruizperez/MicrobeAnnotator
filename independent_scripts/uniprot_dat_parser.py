#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Parses compressed dat files and extracts protein
information relevant for annotation purposes.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import gzip
from pathlib import Path

################################################################################
"""---1.0 Define Functions---"""

def parse_uniprot_dat(dat_file, output_file_table):
    output_folder = Path(output_file_table).parent
    uniprot_to_refseq = output_folder / "03.uniprot_to_refseq.txt"
    with gzip.open(dat_file, 'rt') as uniprot, open(output_file_table, 'w') as output_file, open(uniprot_to_refseq, 'a') as uni_to_ref:
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
                output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene_id, 
                accession, gene_name, ko_number, organism, taxonomy, function, compartment, process, interpro, pfam, ec_number))
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


################################################################################
"""---3.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(description='''This script parses a Uniprot.dat file and output_files a table with\n'''
                                                    '''the ID, Accession, Gene Name, Organism, Taxonomy, KEGG ID, Function,\n
                                                    Compartment, Process, InterPro, and Pfam\n
                                                    For faster usage in alrge files use gnu parallel (read script file to see how)\n'''
                                    '''\nGlobal mandatory parameters: [Input Uniprot.dat File]\n'''
                                    '''\nOptional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='input_dat', action='store', required=True, help='Uniprot.dat file to parse')
    parser.add_argument('-o', '--output_file', dest='output_file_table', action='store', required=False, help='output_file table')
    args = parser.parse_args()

    input_dat = args.input_dat
    output_file_table = args.output_file_table

    parse_uniprot_dat(input_dat, output_file_table)

if __name__ == "__main__":
    main()
