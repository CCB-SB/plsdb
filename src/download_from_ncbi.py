#!/usr/bin/python

import os
import gzip
import pandas
import logging
import argparse
from multiprocessing import Pool

from utils import *

##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    # Input/output
    parser.add_argument('--reports', '-r', help="", required=True, nargs="+")
    parser.add_argument('--seqs', '-e', help="", required=True)
    parser.add_argument('--odir', '-o', help="", required=True)
    parser.add_argument('--suffix', '-s', help="", required=True)
    parser.add_argument('--cores', help='', type=int, default=1)
    return parser

##################################################
# MAIN
##################################################
ODIR = None
REPORTS = None
SEQS = None
SUFFIX = None

def extract_plasmid_features(fname):
    df = read_ncbi_tables([fname])
    # assembly ID
    asm_id = set(df['assembly'])
    assert len(asm_id) == 1, "More than one assembly ID in assembly features file %s: %s" % (fname, ', '.join(asm_id))
    asm_id = asm_id.pop()
    # seq.s
    seq_ids = list((SEQS[SEQS['assembly_accession'] == asm_id])['sequence_accession'])
    # filter
    if len(seq_ids) == 0:
        # no plasmids -> remove
        logging.info("NO PLASMIDS in %s" % fname)
        rm_file(fname)
    else:
        # filter
        df = df[df['genomic_accession'].isin(seq_ids)]
        # save
        compr = None
        if is_gzipped(fname):
            compr = "gzip"
        df.to_csv(path_or_buf=fname, sep="\t", header=True, index=False, index_label=False, compression=compr)

def download_and_extract_plasmids(ftp_path):
    """
    Download assembly file and extract plasmid sequences
    """
    # download
    fname, downloaded = download_ncbi_assembly(ftp_path=ftp_path, suffix=SUFFIX, odir=ODIR, file_can_exist=False)
    # download failed
    if not downloaded:
        rm_file(fname)
        return
    # extract plasmids
    if SUFFIX == 'feature_table.txt.gz':
        try:
            extract_plasmid_features(fname)
        except Exception as e:
            logging.info('EXTR ERROR: %s: %s' % (fname, e))
    return

if __name__ == "__main__":
    # Logger setup
    setup_logger()

    # Args
    args = get_arg_parser().parse_args()

    # Output dir.
    ODIR = args.odir

    # Output file suffix
    SUFFIX = args.suffix

    # Sequence summaries
    SEQS = read_ncbi_tables([args.seqs])
    logging.info("Sequences: %d entries" % SEQS.shape[0])

    # Report files
    REPORTS = read_ncbi_tables(args.reports)
    logging.info("Reports: %d entries" % REPORTS.shape[0])
    # filter
    REPORTS = REPORTS[REPORTS['assembly_accession'].isin(SEQS['assembly_accession'])]
    logging.info("Reports: %d entries were kept" % REPORTS.shape[0])

    # Download
    pool = Pool(args.cores)
    tmp = pool.map(download_and_extract_plasmids, set(REPORTS["ftp_path"]))
    pool.close(); pool.join()
