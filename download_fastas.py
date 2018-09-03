#!/usr/bin/python

import pandas
import logging
import argparse
from multiprocessing import Pool

from utils import split_list, setup_logger, run_cmd

##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tabs', '-t', help="File containing sequence IDs", required=True, nargs="+")
    parser.add_argument('--icol', '-i', help="Column name containing sequence UIDs (same for all files)", required=True)
    parser.add_argument('--ofile', '-o', help="Output file", required=True)
    parser.add_argument('--cores', '-c', help='Number of cores to use', type=int, default=1)
    parser.add_argument('--size', '-s', help='Bulk size', type=int, default=500)
    parser.add_argument('--eutils', '-e', help='Path to dir. with eutils binaries',  default="tools/edirect")
    return parser

##################################################
# MAIN
##################################################
ARGS = None
CMD  = "{epath}/epost -db nuccore -id \"{ids}\" -format uid | {epath}/efetch -format fasta > {ofile}"

def prep_ids(ids):
    """
    [id1, id2, ...] -> "id1,id2,..." with unique IDs
    """
    ids = set(ids)
    return ','.join([str(i)  for i in ids])

def download_sequences(i, ids):
    # tmp output file
    ofile = ARGS.ofile + ".tmp.%d.fasta" % i
    # CMD
    cmd = CMD.format(
        epath=ARGS.eutils,
        ids=prep_ids(ids),
        ofile=ofile
    )
    logging.info('CMD: %s' % cmd)
    cmd, cmd_s, cmd_o = run_cmd(cmd)
    assert cmd_s == 0, "CMD ERROR: %s: %d\n%s" % (cmd, cmd_s, cmd_o)
    return

if __name__ == "__main__":
    # Logger setup
    setup_logger()

    # Args
    ARGS = get_arg_parser().parse_args()

    # Tables
    ids = set()
    for tab in ARGS.tabs:
        tmp = pandas.read_csv(filepath_or_buffer=tab, sep='\t', header=0)
        assert ARGS.icol in list(tmp.columns), 'File %s has no column \"%s\"' % (tab, ARGS.icol)
        logging.info('TAB %s: %d IDs' % (tab, tmp.shape[0]))
        ids.update(list(tmp[ARGS.icol]))
    logging.info('IDS: %d' % len(ids))

    # Download
    pool = Pool(ARGS.cores)
    tmp = pool.starmap(
        download_sequences,
        enumerate([chunk for chunk in split_list(ids, ARGS.size)])
    )
    pool.close()
    pool.join()
