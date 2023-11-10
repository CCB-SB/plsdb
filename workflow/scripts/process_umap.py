import argparse
import logging
import pandas
import umap
import numpy
from utils_my_logger import My_logger



##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help="Input", required=True)
    parser.add_argument('--ofile', '-o', help="Output file.", required=True)
    parser.add_argument('--log', '-l', help="Log file.", required=True)
    parser.add_argument('--neighbors', help="UMAP neighbors.", type= int, required=True)
    parser.add_argument('--components', help="UMAP components.", type=int, required=True)
    parser.add_argument('--min_dist', help="UMAP minimal distance.", type=float, required=True)
    return parser


##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    ARGS = get_arg_parser().parse_args()
    logger = My_logger(log_filename = ARGS.log, logger_name = "process_umap")
    logger = logger.get_logger()

    # IDs
    ids = None
    with open(ARGS.input, 'r') as ifile:
        ids = ifile.readline().rstrip('\n').split('\t')[1:]
        logger.info('Distance matrix: {num} x {num}'.format(num=len(ids)))

    # distance matrix
    logger.info('Start loading distances...')
    dist = numpy.loadtxt(
        fname=ARGS.input,
        skiprows=1,
        usecols=range(1,len(ids) + 1) # skip 1st (contains IDs)
    )

    # embedding
    logger.info('Start embedding...')
    embedding = umap.UMAP(
        n_neighbors=ARGS.neighbors,
        n_components=ARGS.components,
        min_dist=ARGS.min_dist,
        init='random',
        metric='precomputed',
        random_state=42
    ).fit_transform(dist)
    logger.info('Done.')

    # save to file
    embedding = pandas.DataFrame(
        embedding,
        columns=['D1', 'D2'],
        index=ids
    )
    embedding.to_csv(
        path_or_buf=ARGS.ofile,
        sep='\t',
        header=True,
        index=True,
        index_label='ID'
    )