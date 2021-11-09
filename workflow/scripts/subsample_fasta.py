from Bio import SeqIO
import argparse
import logging
import pandas
import os
from tqdm import tqdm
from utils import setup_logger, mkdir

##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ifile', '-i', help="Input tsv file with ACC_NUCCORE column", required=True)
    parser.add_argument('--fasta', '-f', help="File containing sequence IDs", required=True)
    parser.add_argument('--ofile', '-o', help="Output file", required=True)
    parser.add_argument('--log', '-l', help="logfile", required=False)
    return parser





##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    # Args
    ARGS = get_arg_parser().parse_args()

    logger = setup_logger(logging.INFO,ARGS.log)

    mkdir(os.path.dirname(ARGS.ofile), True)

    ids = list(pandas.read_csv(ARGS.ifile, sep='\t', header=0, dtype=str)['ACC_NUCCORE'])
    n   = len(ids)
    ids = set(ids)
    assert len(ids) == n, 'FASTA IDs in {} are not unique'.format(ARGS.ifile)
    logger.info('Read in {} plasmid records'.format(len(ids)))

    with open(ARGS.fasta, 'r') as ifile, open(ARGS.ofile, 'w') as ofile:
        for record in tqdm(SeqIO.parse(ifile, 'fasta')):
            if record.id in ids:
                SeqIO.write(record, ofile, 'fasta')
                ids.remove(record.id)
    assert len(ids) == 0, 'Not all IDs found in FASTA {}: {}'.format(ARGS.fasta, ';'.join(ids))