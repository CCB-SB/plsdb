import argparse
import pandas
from utils_my_logger import My_logger
from utils_decorators import timer, debuger

import logging
logg = logging.getLogger("rmlst_blastn")
##################################################
# FUNCTIONS
##################################################


##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', '-f', help="Fasta file.", required=True)
    parser.add_argument('--ofile', '-o', help="Output file.", required=True)
    parser.add_argument('--db', help="rMLST fasta.", required=True)
    parser.add_argument('--dbs', help="rMLST db.", required=True, nargs='+')
    parser.add_argument('--cores', help="Number of cores for this rule.", required=True)
    parser.add_argument('--header', help="rMLST header from config file.", required=True, nargs='+')
    parser.add_argument('--ident', help="rMLST ident from config file.", required=True)
    parser.add_argument('--cov', help="rMLST cov.", type=float, required=True)
    parser.add_argument('--log', help="Log file.", required=True)
    return parser


##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    from utils import run_cmd
    ARGS = get_arg_parser().parse_args()

    # Logger
    logger = My_logger(log_filename = ARGS.log, logger_name = "rmlst_blastn")
    logger = logger.get_logger()

    # blastn search
    logger.info(f"Searching plasmid sequences against rmlst database using Blastn...")
    cmd = f"""blastn -query {ARGS.fasta} -db {ARGS.db} -task blastn -perc_identity {ARGS.ident}\
        -out {ARGS.ofile} -outfmt '6 {' '.join(ARGS.header)}' -num_threads {ARGS.cores}\
        -evalue 1E-20 -culling_limit 1"""
    logger.info(cmd)
    run_cmd(cmd, split=True, shell=False)
    logger.info("Blastn completed")

    # read in results
    df = pandas.read_csv(ARGS.ofile, sep='\t', header=None, names=ARGS.header)
    logger.info("Blastn output file read")

    # compute coverage
    df['cov'] =  100 * (df['length'] - df['gaps']) / df['slen']
    
    # filter by coverage
    df = df.loc[df['cov'] >= ARGS.cov,:]

    # rMLST locus
    df['slocus'] = df['sseqid'].map(lambda x: x.split('_')[0])

    # save
    df.to_csv(
        path_or_buf=ARGS.ofile,
        sep='\t',
        header=True,
        index=False,
        index_label=False
    )

    logger.info(f"File generated: {ARGS.ofile}")
