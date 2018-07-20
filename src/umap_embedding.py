#!/usr/bin/python

import os
import umap
import numpy
import pickle
import pandas
import logging
import argparse

from utils import *

##################################################
# GLOBAL
##################################################
# XXX

##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    # Input/output
    parser.add_argument('--dist', '-d', help="Distance matrix (Mash output)", required=True)
    parser.add_argument('--ofile', '-o', help="Output file", required=True)
    # UMAP params
    parser.add_argument('--n_neighbors', default=50, type=int)
    parser.add_argument('--n_components', default=2, type=int)
    return parser

##################################################
# MAIN
###############################################
if __name__ == "__main__":
    # Logger setup
    setup_logger()

    # Args
    args = get_arg_parser().parse_args()

    # Distances
    ncols = 0
    with open(args.dist) as f:
        ncols = len(f.readline().split('\t'))
        logging.info('There are %d columns' % ncols)
    dist = pandas.read_csv(
        args.dist, sep='\t',
        header=None, index_col=None,
        dtype=numpy.float16,
        skiprows=1, usecols=range(1, ncols + 1),
        low_memory=True
    )
    logging.info('Distance matrix: %d x %d' % (dist.shape[0], dist.shape[1]))

    # Sample names
    ids  = dist.columns

    # UMAP
    logging.info('Start embedding')
    embedding = umap.UMAP(
        n_neighbors=args.n_neighbors,
        n_components=args.n_components,
        init='random',
        metric='precomputed'
    ).fit_transform(dist.as_matrix())
    logging.info('Finished')

    # Convert to DataFrame
    embedding = pandas.DataFrame(
        embedding,
        columns=['D1', 'D2'],
        index=dist.columns
    )

    # Write to file
    logging.info('Save to %s' % args.ofile)
    embedding.to_csv(
        path_or_buf=args.ofile,
        sep='\t',
        header=True,
        index=True,
        index_label='ID'
    )
