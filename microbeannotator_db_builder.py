#!/usr/bin/env python

"""
########################################################################
# Author:       Carlos A. Ruiz Perez
# Email:        cruizperez3@gatech.edu
# Intitution:   Georgia Institute of Technology
# Version:      1.0
# Date:         Nov 13, 2020

# Description: Downloads all data required to build the search databases.
# Parses the annotation information and creates the method-specific dbs.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
from independent_scripts import conversion_database_creator
from independent_scripts import uniprot_data_downloader
from independent_scripts import refseq_data_downloader
from independent_scripts import uniprot_dat_parser
from independent_scripts import refseq_genbank_parser
from independent_scripts import sqlite3_database_builder
from independent_scripts import protein_db_creation
from pathlib import Path
from shutil import rmtree


################################################################################
"""---1.0 Define Functions---"""
def db_duilder(output_file_folder, method, light, threads, step, aspera, keep, path=None):
    # Create output_file folder
    outfolder = Path(output_file_folder)
    outfolder.mkdir(parents=True, exist_ok=True)
    genbank_files = None
    database = outfolder / "02.MicrobeAnnotator.db"
    # Download information
    if step == 1:
        print("Downloading UniProt information...")
        uniprot_data_downloader.uniprot_fasta_downloader(outfolder, light)
        uniprot_data_downloader.uniprot_dat_downloader(outfolder, light)
        print("Finished downloading UniProt data!")
        if light == False:
            if aspera == True:
                print("Downloading RefSeq information using Aspera Connect...")
                refseq_data_downloader.refseq_fasta_downloader(outfolder)
                genbank_files = refseq_data_downloader.refseq_genbank_downloader(outfolder)
            else:
                print("Downloading RefSeq information using wget...")
                refseq_data_downloader.refseq_fasta_downloader_wget(outfolder, threads)
                genbank_files = refseq_data_downloader.refseq_genbank_downloader_wget(outfolder, threads)
            print("Finished downloading RefSeq data!")
        step += 1
    # -------------------

    # Parse annotation files
    if step == 2:
        if light == True:
            print("Parsing swissprot annotation data...")
            temporal_sprot_dat = outfolder / "02.temp_dat_files/uniprot_sprot.dat.gz"
            final_sprot_table = outfolder / "02.uniprot_sprot.table"
            uniprot_dat_parser.parse_uniprot_dat(temporal_sprot_dat, final_sprot_table)
            # Remove temporary files
            if keep == False:
                if Path(outfolder / "02.temp_dat_files").is_dir():
                    rmtree(outfolder / "02.temp_dat_files")
            print("Done!")
        else:
            if genbank_files == None and aspera == True:
                temp_gb_files = Path(output_file_folder) / "02.temp_genbank"
                viral_list = refseq_data_downloader.search_all_files(temp_gb_files / "viral")
                bacteria_list = refseq_data_downloader.search_all_files(temp_gb_files / "bacteria")
                archaea_list = refseq_data_downloader.search_all_files(temp_gb_files / "archaea")
                genbank_files = viral_list + bacteria_list + archaea_list
            elif genbank_files == None and aspera == False:
                temp_gb_files = Path(output_file_folder) / "02.temp_genbank"
                genbank_files = refseq_data_downloader.search_all_files(temp_gb_files)
            print("Parsing swissprot and trembl annotation data...")
            temporal_sprot_dat = outfolder / "02.temp_dat_files/uniprot_sprot.dat.gz"
            final_sprot_table = outfolder / "02.uniprot_sprot.table"
            uniprot_dat_parser.parse_uniprot_dat(temporal_sprot_dat, final_sprot_table)
            temporal_trembl_dat = outfolder / "02.temp_dat_files/uniprot_trembl.dat.gz"
            final_trembl_table = outfolder / "02.uniprot_trembl.table"
            uniprot_dat_parser.parse_uniprot_dat(temporal_trembl_dat, final_trembl_table)

            print("Done!")
            print("Parsing refseq annotation data...")
            temp_table_list = refseq_genbank_parser.table_generator_worker(genbank_files, threads)
            final_refseq_table = outfolder / "02.refseq.table"
            refseq_genbank_parser.table_merger(temp_table_list, final_refseq_table, False)
            # Remove temporary files
            if keep == False:
                if Path(outfolder / "02.refseq.table").is_dir():
                    rmtree(outfolder / "02.refseq.table")
                if Path(outfolder / "02.temp_dat_files").is_dir():
                    rmtree(outfolder / "02.temp_dat_files")
            print("Done!")
        step += 1
    # --------------------

    # Merge annotations into SQLite database
    if step == 3:
        print("Building SQLite database: 02.MicrobeAnnotator.db")
        if light == True:
            final_sprot_table = outfolder / "02.uniprot_sprot.table"
            sqlite3_database_builder.create_swissprot_table(str(database), final_sprot_table)
            if keep == False:
                Path(final_sprot_table).unlink()
        else:
            final_sprot_table = outfolder / "02.uniprot_sprot.table"
            final_trembl_table = outfolder / "02.uniprot_trembl.table"
            final_refseq_table = outfolder / "02.refseq.table"
            sqlite3_database_builder.create_swissprot_table(str(database), final_sprot_table)
            sqlite3_database_builder.create_trembl_table(str(database), final_trembl_table)
            sqlite3_database_builder.create_refseq_table(str(database), final_refseq_table)
            if keep == False:
                if Path(final_sprot_table).is_file():
                    Path(final_sprot_table).unlink()
                if Path(final_trembl_table).is_file():
                    Path(final_trembl_table).unlink()
                if Path(final_refseq_table).is_file():
                    Path(final_refseq_table).unlink()
        print("Done!")
        step += 1
    # ---------------------

    # Create interconversion tables
    if step == 4:
        print("Building interconversion databases...")
        interconversion_database = outfolder / "03.Conversion.db"
        print("Building RefSeq to Uniprot")
        try:
            refseq_to_uniprot_table = str(outfolder / "03.uniprot_to_refseq.txt")
            conversion_database_creator.create_refseq_to_uniprot(refseq_to_uniprot_table, interconversion_database, keep)
        except:
            print("Could not create RefSeq to Uniprot table.")
        print("Building KO to EC")
        try:
            conversion_database_creator.create_ko_to_ec(outfolder, interconversion_database, keep)
        except:
            print("Could not create KO to EC table.")
        print("Building Interpro to EC and KO")
        try:
            conversion_database_creator.create_interpro_tables(outfolder, interconversion_database, keep)
        except:
            print("Could not create UniProt tables.")
        step += 1
        print("Done!")
    
    # Create method-specific protein databases
    if step == 5:
        print("Building protein databases...")
        protein_folder = outfolder / "01.Protein_DB"
        if method == "blast":
            protein_db_creation.blastp_db_creator(protein_folder, path)
        elif method == "diamond":
            protein_db_creation.diamond_db_creator(protein_folder, threads, path)
        elif method == "sword":
            protein_db_creation.sword_db_creator(protein_folder, path)
        else:
            print("Please provide a valid db creation method.")
            print("If you are running the microbeannotator_db_builder script this is the last step of the process.")
            print('Given that this last step failed you need to re-run this script starting from --step 4')
            print('and provide the binary folder for the selected method using the "--bin_path" option or add')
            exit('the binary to your PATH (see README for help on this).')
        print("Done!")
        print("---------")
        print("MicrobeAnnotator has finished creating your databases! You are ready to search and annotate your proteins!")


################################################################################
"""---2.0 Main Function---"""

def main():
    import argparse, sys
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script build the search databases required by MicrobeAnnotator\n'''
            '''Usage: ''' + sys.argv[0] + ''' -f [output_file folder]\n'''
            '''Global mandatory parameters: -f [output_file folder]\n'''
            '''Optional Database Parameters: See ''' + sys.argv[0] + ' -h')
    parser.add_argument('-d', '--dir', dest='directory', action='store', required=True,
                        help='Directory where database will be created.')
    parser.add_argument('-m', '--method', dest='method', action='store', required=True,
                        help='Search (and db creation) method, one of blast, diamond or sword')
    parser.add_argument('-t', '--threads', dest='threads', action='store', required=False, default=1, type=int,
                        help='Threads to use (when possible). By default 1.')
    parser.add_argument('--bin_path', dest='bin_path', action='store', required=False,
                        help='Path to binary folder for selected method. By defaul assumes the program is in path.')
    parser.add_argument('--step', dest='step', action='store', required=False, default=1, type=int,
                        help='Step to start with. 1.Download data, 2.Parse annotation data, 3.Building SQLite DB, 4.Create interconversion tables, 5.Build protein DBs. Default 1.')
    parser.add_argument('--light', dest='light', action='store_true', required=False,
                        help='Use only KOfamscan and swissprot databases. By default also builds refseq and trembl.')
    parser.add_argument('--no_aspera', dest='aspera', action='store_false', required=False,
                        help='Disables download using Aspera and instead uses wget. By default uses Aspera Connect.')
    parser.add_argument('--keep', dest='keep', action='store_true', required=False,
                        help='Keep intermediate files, can increase disk requirement (not necessary and therefore not recommended)')
    args = parser.parse_args()

    directory = args.directory
    method = args.method
    method = method.lower()
    bin_path = args.bin_path
    light = args.light
    threads = args.threads
    step = args.step
    aspera = args.aspera
    keep = args.keep

    # ----------------------------
    db_duilder(directory, method, light, threads, step, aspera, keep, bin_path)
    # ----------------------------

if __name__ == "__main__":
    main()
