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
import re
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd

##################################################
# FUNCTIONS
##################################################




##################################################
# ARGS
##################################################

def get_parser():
    import argparse
    import json
    # create parser object
    parser = argparse.ArgumentParser(
        description="Assign BIOSAMPLE Host using iso/host mappings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # specify our desired options 
    # by default ArgumentParser will add an help option 
    parser.add_argument("--input-pls", help="csv file for plasmid metadata")
    parser.add_argument("--input-fasta", help="Fasta file of plasmid metadata")
    parser.add_argument("--threads", type = int, default = 1, required = False,
        help="Threads to use" )
    parser.add_argument("--log")
    parser.add_argument("--outfile")

    return parser

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    ARGS = get_parser().parse_args()
    
    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "process_calculate_gc")
    my_logger = logger.get_logger()

    #
    ## Here starts
    #

    # Read plasmid md
    pls = pd.read_csv(ARGS.input_pls, header=0, dtype=str)
    pls.set_index('NUCCORE_ACC', drop=False, inplace=True, verify_integrity=True)
    my_logger.info('Read in {} plasmid records'.format(pls.shape[0]))

    # collect seq. stats.s
    my_logger.info('Collecting sequences and their data ...')

    # Read fasta
    seq = []
    seqs = {}
    with open(ARGS.input_fasta, 'r') as ifile:
        for record in SeqIO.parse(ifile, 'fasta'):
            assert re.fullmatch(r'\w+\.\d+', record.id), 'Unexpected ID format in FASTA: {}'.format(record.id)
            assert record.id in list(pls.index), 'Unknown FASTA ID {}'.format(record.id)
            seq.append({
                'NUCCORE_ACC': record.id,
                'NUCCORE_GC': gc_fraction(record.seq),
                'Length': len(record.seq)
            })
            seqs[record.id] = record.seq
    
    seq = pd.DataFrame(seq)
    seq.set_index('NUCCORE_ACC', drop=True, inplace=True, verify_integrity=True)

    # add to table
    pls = pd.merge(
        left=pls,
        right=seq,
        how='left',
        left_index=True,
        right_index=True,
        sort=False,
    )

    # check seq. length
    if pls['NUCCORE_Length'].equals(pls['Length']):
        my_logger.warning('Sequence length is not always identical: nuccore vs. FASTA')

    # save
    pls.to_csv(ARGS.outfile, header=True, index=False, index_label=False)
