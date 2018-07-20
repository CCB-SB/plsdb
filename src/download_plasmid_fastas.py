#!/usr/bin/python

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
    parser.add_argument('--odir', '-o', help="", required=True)
    parser.add_argument('--suffix', '-s', help="", default="genomic.fna.gz")
    parser.add_argument('--cores', help='', type=int, default=1)
    return parser

##################################################
# MAIN
##################################################
ODIR = None
REPORTS = None
SUFFIX = None

def extract_plasmid_seqs(fname):
    """
    Filter given FASTA for plasmid sequences
    """
    plasmids = []
    for i, record in enumerate(get_fasta_records(fname)):
        record_type = get_record_type(record)
        logging.info('TYPE: %s: %s in %s' % (record_type, record.id, fname))
        if record_type == 'plasmid':
            logging.info('PLASMID: %s in %s' % (record.description, fname))
            plasmids.append(record)
    if len(plasmids) == 0:
        # no plasmids -> remove
        logging.info("NO PLASMIDS in %s" % fname)
        rm_file(fname)
    else:
        records2fasta(plasmids, fname)

def download_and_extract_plasmids(ftp_path):
    """
    Download FASTA file and extract plasmid sequences
    """
    # download
    fname, downloaded = download_ncbi_assembly(ftp_path=ftp_path, suffix=SUFFIX, odir=ODIR, file_can_exist=False)
    # download failed
    if not downloaded: # will log errors
        rm_file(fname)
        return
    # extract plasmids
    try:
        extract_plasmid_seqs(fname)
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

    # Report files
    REPORTS = pandas.concat([
        pandas.read_csv(filepath_or_buffer=report_file, sep='\t', header=0)
        for report_file in args.reports
    ])
    logging.info("Reports: %d entries" % REPORTS.shape[0])

    # Download
    pool = Pool(args.cores)
    tmp = pool.map(download_and_extract_plasmids, set(REPORTS["ftp_path"]))
    pool.close(); pool.join()
