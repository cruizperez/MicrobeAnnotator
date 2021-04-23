#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos A. Ruiz-Perez
# Email:       cruizperez3@gatech.edu
# Institution: Georgia Institute of Technology
# Version:	   0.9
# Date:		   14 February 2020

# Description: This script filters a fasta file based on a list of IDs
# to retain or exclude
########################################################################
"""

################################################################################
"""---1.0 Import Modules---"""
import argparse
from sys import argv

################################################################################
"""---2.0 Define Functions---"""

def fastA_filter_list(fasta_file, output_file, id_list, reverse=False):
    # import modules
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    # Filter fasta file
    with open(fasta_file) as fasta_input, open(output_file, 'w') as output:
        if reverse == True:
            for title, seq in SimpleFastaParser(fasta_input):
                if title.split()[0] not in id_list:
                    output.write(">{}\n{}\n".format(title, seq))
        else:
            for title, seq in SimpleFastaParser(fasta_input):
                if len(id_list) < 1:
                    break
                elif title.split()[0] in id_list:
                    output.write(">{}\n{}\n".format(title, seq))
                    id_list.remove(title.split()[0])

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
    description='''Filter a FastA file based ona provided list of IDs (file or input).\n'''
    '''It can exclude or retrieve the sequences using the --reverse flag\n'''
    '''Usage: ''' + argv[0] + ''' -f [FastA File] -o [output_file File] -l [File with IDs] OR -i [ID list]\n'''
    '''Global mandatory parameters: [FastA_File] [output_file_File] [ID List File]\n'''
    'Optional Database Parameters: See ' + argv[0] + ' -h')
    parser.add_argument('-f', '--fasta', dest='fasta_file', action='store', required=True, 
    help='FastA file to filter')
    parser.add_argument('-o', '--output_file', dest='output_file', action='store', required=True, 
    help='output_file FastA file with retrieved sequences')
    parser.add_argument('-l', '--list', dest='id_file', action='store', required=False, 
    help='File with list of IDs to filter.')
    parser.add_argument('-i', '--id', dest='id_list', action='store', required=False,
    help='Comma-separated IDs to filter: ID1,ID2,ID3')
    parser.add_argument('--reverse', action='store_true', 
    help='Exclude the sequences in the list file. By default False, i.e. retrieves those in the list')
    args = parser.parse_args()

    fasta_file = args.fasta_file
    output_file = args.output_file
    id_file = args.id_file
    id_list = args.id_list
    reverse = args.reverse

    if id_list != None:
        id_list = id_list.split(",")
        fastA_filter_list(fasta_file, output_file, id_list, reverse)
    elif id_file != None:
        id_list = []
        with open(id_file, 'r') as file:
            for line in file:
                id_list.append(line.strip())
        fastA_filter_list(fasta_file, output_file, id_list, reverse)
    else:
        raise ValueError("Did you provide the IDs to filter?")

if __name__ == "__main__":
    main()