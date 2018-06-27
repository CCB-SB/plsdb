#!/usr/bin/python

import os
import pandas
import logging
import argparse
from glob import glob

from utils import *

##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    # Input/output
    parser.add_argument('--fdir', '-f', help="", required=True)
    parser.add_argument('--ofile', '-o', help="", required=True)
    return parser

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    # Logger setup
    setup_logger()

    # Args
    args = get_arg_parser().parse_args()

    # Suffix
    suffix = "genomic.fna.gz"

    # Downloaded files
    dfiles = sorted(glob("%s/*_%s" % (args.fdir, suffix)))
    logging.info("Found %d files in %s with *_%s" % (len(dfiles), args.fdir, suffix))

    # Suffix for re
    re_suffix = re.compile("_" + re.sub("\.", "\.", suffix))

    # record ID -> assembly ID
    rinfo = {}
    dupl  = {}
    for dfile in dfiles:
        asm_id = os.path.basename(dfile)
        asm_id = re.sub(re_suffix, "", asm_id)
        asm_id = (
            '_'.join(asm_id.split('_')[0:2]),
            '_'.join(asm_id.split('_')[2:])
        )
        # for each sequence in FASTA
        for record in fasta_records(dfile):
            # already saved
            if record.id in rinfo:
                dupl[record.id] = dfile
            rinfo[record.id] = {
                'assembly_accession': asm_id[0],
                'asm_name': asm_id[1],
                'sequence_description': ' '.join(record.description.split(' ')[1:])
            }
    logging.info("Found %d records in %d files" % (len(rinfo), len(dfiles)))
    assert len(dupl) == 0, "There are %d duplicate records" % len(dupl)

    # Record table
    rinfo = pandas.DataFrame.from_dict(rinfo, orient="index")
    # add row names as column
    rinfo = rinfo.assign(sequence_accession=rinfo.index)
    # re-order columns
    rinfo = rinfo[["sequence_accession", "sequence_description", "assembly_accession", "asm_name"]]

    # save
    rinfo.to_csv(path_or_buf=args.ofile, sep="\t", header=True, index=False, index_label=False)
