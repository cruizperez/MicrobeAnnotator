#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0.0
# Date:         Nov 13, 2020

# Description: Creates a table to KO to EC table based on KEGG information.
# This script requires you to download the 'htext' file from
# https://www.kegg.jp/kegg-bin/get_htext?ko00001
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""


################################################################################
"""---1.0 Define Functions---"""
def kegg_orthology_parser(kegg_file, outfile):
    ko_ec_dict= {}
    with open(kegg_file, 'r') as infile:
        for line in infile:
            if line.startswith('D') and "EC:" in line:
                ko_id = line.strip().split()[1]
                ec_numbers = line.strip().split("EC:")[1]
                ec_numbers = ec_numbers.replace("]", "")
                ec_numbers = ec_numbers.split()
                ko_ec_dict[ko_id] = ec_numbers
    with open(outfile, 'w') as output:
        for ko_identifier, ec_identifiers in ko_ec_dict.items():
            for element in ec_identifiers:
                output.write("{}\t{}\n".format(ko_identifier, element))

################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''Creates a table to KO to EC table based on KEGG information.\n'''
            '''This script requires you to download the 'htext' file from\n'''
            '''https://www.kegg.jp/kegg-bin/get_htext?ko00001.\n'''
            '''Usage: ''' + sys.argv[0] + ''' -i [input_table] -o [output_table]\n'''
            '''Global mandatory parameters: -i [input_table] -o [output_table]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-i', '--input', dest='infile', action='store', required=False,
                        help='Input table downloaded from KEGG')
    parser.add_argument('-o', '--output', dest='output', action='store', required=False,
                        help='Output table with ko to ec ids.')
    args = parser.parse_args()

    infile = args.infile
    output = args.output

    # ----------------------------
    kegg_orthology_parser(infile, output)
    # ----------------------------

if __name__ == "__main__":
    main()
