#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 
# Authors:   Alejandra Gonzalez (gonmola@hotmail.es)
# ---


##################################################
# IMPORTS
##################################################


from utils_my_logger import My_logger
from utils_decorators import timer, debuger

import logging
logg = logging.getLogger("plasmid_sequence_filtering")


##################################################
# FUNCTIONS
##################################################

def remove_duplicates(input_file, outfile, my_logger=None):
    import subprocess
    my_logger.info(f"Checking for possible duplicated sequences in {input_file}")

    tmp_file = "tmp.duplicated_ids.txt" # info of duplicated IDs
    cmd = f"cat {input_file} | seqkit rmdup --by-seq --dup-num-file {tmp_file} > {outfile}"
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    proc.wait()
    message = proc.stderr.read().decode()

    my_logger.info(f"Seqkit output: {message}")

    return tmp_file

##################################################
# ARGS
##################################################
def get_arg_parser():
    import argparse

    parser = argparse.ArgumentParser(description="Filter duplicated sequences",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Mandatory argument
    parser.add_argument("--fasta-file", required = True, 
                        help="Fasta file")
    parser.add_argument("--metadata-file", required = True, 
                        help="CSV file")
    # Optional args
    parser.add_argument("--outfile-seqs", default = "seq_filtered.fna", required = False,
                        help="Output file for sequences" )
    parser.add_argument("--outfile-md", default = "md_filtered.csv", required = False,
                        help="Output file for metadata CSV" )
    parser.add_argument("--log", default = "log.log", required = False)

    return parser

##################################################
# MAIN
##################################################

if __name__ == "__main__":
    
    from os.path import exists, join
    import os
    import shutil
    import numpy as np
    import pandas as pd
    from utils import run_cmd, filter_fasta
    ARGS = get_arg_parser().parse_args()
    
    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "plasmid_sequence_filtering")
    my_logger = logger.get_logger()

    # Look for duplicated sequences
    tmp_outfile = "tmp.sequences.fna"
    tmp_file = remove_duplicates(input_file = ARGS.fasta_file, outfile = tmp_outfile, my_logger=my_logger)
    if not exists(tmp_file):
        my_logger.info("NO Duplicated sequences detected.")
        os.rename(tmp_outfile, ARGS.outfile_seqs)
        shutil.copy2(ARGS.metadata_file, ARGS.outfile_md)
        my_logger.info("Generated metadata and sequence files are copies from input files")
        exit()
    else:
        my_logger.info("Duplicated sequences detected.")

    # Load metadata file
    md = pd.read_csv(ARGS.metadata_file)
    ## Date format
    date_cols = ["ASSEMBLY_SeqReleaseDate","NUCCORE_CreateDate" ]
    for col in date_cols:
        md[col] = pd.to_datetime(
            md[col],
            format='%Y-%m-%d')
    ## Check if Assembly and Biosample Info
    md['has_biosample'] = np.where(md['BIOSAMPLE_ACC'].isnull(), False, True)
    md['has_assembly'] = np.where(md['ASSEMBLY_ACC'].isnull(), False, True)
    md['has_location'] = np.where(md['BIOSAMPLE_Location'].isnull(), False, True)

    # Parse duplicated seqs id info in list of tuples
    dup_groups = [] 
    with open(tmp_file, 'r') as f:
        for line in f:
            # [NUCCORE_ACC, NUCCORE_ACC,...,]
            nuccore_acc_list = line.strip().split("\t")[1].split(', ')
            dup_groups.append(nuccore_acc_list)
    
    
    selecting_criteria = """
    Compare two plasmid records:
    1) Prefer the one from RefSeq
    2) Prefer the one with location information
    3) Prefer the one with an assembly
    4) Prefer the most recent assembly release date
    5) Prefer the one with a biosample
    6) Prefer the newewst nuccore creation date
    7) Prefer the one with the highest coverage
    8) If all equals, choose the first one
    """
    sorting_columns = [
            'NUCCORE_Source', "has_location","has_assembly", 
            "ASSEMBLY_SeqReleaseDate", "has_biosample",
            "NUCCORE_CreateDate", "ASSEMBLY_coverage"
            ]
    my_logger.info(f"""There are {len(dup_groups)} groups of equal sequences.\nTo choose the nuccore record,
                   the criteria will be the following:\n {selecting_criteria}""")

    # Compare records
    nequals = 0
    ids_to_delete = []
    for group in dup_groups:
        # subset md
        md_sub = md[md['NUCCORE_ACC'].isin(group)]
        
        # Sort by criteria
        md_sub = md_sub.sort_values(by=sorting_columns,
            ascending=[False, False, False, False, False, False, False])

        # Check if first and second items are equals
        record_to_keep = md_sub.iloc[0]
        equals = record_to_keep[sorting_columns].equals(md_sub[sorting_columns].iloc[1])

        if equals:
            nequals = nequals + 1
        
        # Select the the first item and keep record of the others
        records_to_delete = md_sub[md_sub['NUCCORE_ACC']!= record_to_keep['NUCCORE_ACC']]
        to_delete = list(records_to_delete['NUCCORE_ACC'])
        ids_to_delete.extend(to_delete)
    
    my_logger.info(f"Sequences equals although criteria: {nequals}")
    
    # Filter fasta
    nrecords = filter_fasta(
        ID_list = ids_to_delete, 
        input_fasta = ARGS.fasta_file, out_fasta = ARGS.outfile_seqs,
        tmp_file = tmp_file, grep_invert = True)
    my_logger.info(f"Filtered fasta file: {nrecords} records  - {ARGS.outfile_seqs} ")
    
    # Filter metadata
    md = md[~md['NUCCORE_ACC'].isin(ids_to_delete)]
    md.to_csv(ARGS.outfile_md, index=False)
    my_logger.info(f"Filtered metadata file: {len(md.index)} records - {ARGS.outfile_md}")
    
    # Remove temporal files
    os.remove(tmp_outfile)


