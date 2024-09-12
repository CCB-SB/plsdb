#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 
# Authors:   Alejandra Gonzalez (gonmola@hotmail.es)
# ---

##################################################
# IMPORTS
##################################################
from utilsmeta.my_logger import My_logger
from utilsmeta.utils import load_table
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import re
import pandas as pd

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    
    # Logging
    logger = My_logger(log_filename = snakemake.log[0], logger_name = "process_calculate_gc")
    my_logger = logger.get_logger()

    #
    ## Here starts
    #

    # Read plasmid md
    pls = load_table(snakemake.input.pls, table="NUCCORE")
    pls.set_index('NUCCORE_ACC', drop=False, inplace=True, verify_integrity=True)
    my_logger.info(f'Read in {len(pls.index)} plasmid records')

    # collect seq. stats.s
    my_logger.info('Collecting sequences and their data ...')

    # Read fasta
    seq = []
    seqs = {}
    with open(snakemake.input.fasta, 'r') as ifile:
        for record in SeqIO.parse(ifile, 'fasta'):
            assert re.fullmatch(r'\w+\.\d+', record.id), f'Unexpected ID format in FASTA: {record.id}'
            assert record.id in list(pls.index), f'Unknown FASTA ID {record.id}'
            seq.append({
                'NUCCORE_ACC': record.id,
                'NUCCORE_GC': gc_fraction(record.seq),
                'Length': len(record.seq)
            })
            seqs[record.id] = record.seq
    
    seq = pd.DataFrame(seq)
    seq.set_index('NUCCORE_ACC', drop=True, inplace=True, verify_integrity=True)

    # add to table
    pls = pd.merge(left=pls, right=seq, how='left',
                   left_index=True, right_index=True, sort=False)

    # check seq. length
    if not pls['NUCCORE_Length'].equals(pls['Length']):
        my_logger.warning('Sequence length is not always identical: nuccore vs. FASTA')

    # save
    pls.to_csv(snakemake.output[0], header=True, index=False, index_label=False)
