#!/usr/bin/python

import pandas
import argparse
from Bio import SeqIO
from multiprocessing import Pool
import time

from utils import split_list, run_cmd

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
    parser.add_argument('--eutils', '-e', help='Path to dir. with eutils binaries', default="tools/edirect")
    parser.add_argument('--idformat', '-f', help='ID format: uid or acc', default="uid", choices=['uid', 'acc'])
    return parser

##################################################
# MAIN
##################################################
ARGS = None
CMD  = "{epath}/epost -db nuccore -id \"{ids}\" -format {IDformat} | {epath}/efetch -format fasta > {ofile}"

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
        ofile=ofile,
        IDformat=ARGS.idformat
    )
    # Execute
    
    count = 0
    for repeat in range (5):
        cmd, cmd_s, cmd_o = run_cmd(cmd)
        # Check CMD status
        assert cmd_s == 0, "CMD ERROR: %s: %d\n%s" % (cmd, cmd_s, cmd_o)
        # Check FASTA
        with open(ofile, 'r') as of:
            count = 0
            for record in SeqIO.parse(of, 'fasta'):
                count += 1
            if (count == len(set(ids))): # length is correct
                return
            # else: repeat the loop and therefore: rerun!
            time.sleep(5) # sleep in case NCBI has problems with too many accesses
    assert count == len(set(ids)), "FASTA ERROR: missing IDs: %s" % (cmd)
    return

if __name__ == "__main__":
    # Args
    ARGS = get_arg_parser().parse_args()

    # Tables
    ids = set()
    for tab in ARGS.tabs:
        tmp = pandas.read_csv(filepath_or_buffer=tab, sep='\t', header=0)
        assert ARGS.icol in list(tmp.columns), 'File %s has no column \"%s\"' % (tab, ARGS.icol)
        print('TAB %s: %d IDs' % (tab, tmp.shape[0]))
        ids.update(list(tmp[ARGS.icol]))
    print('IDS: %d' % len(ids))

    # Download
    pool = Pool(ARGS.cores)
    tmp = pool.starmap(
        download_sequences,
        enumerate([chunk for chunk in split_list(ids, ARGS.size)])
    )
    pool.close()
    pool.join()
    