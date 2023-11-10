#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 
# Authors:   Alejandra Gonzalez (gonmola@hotmail.es)
# ---

###########
# IMPORTS #
###########

from utils_my_logger import My_logger
from utils_decorators import timer, debuger
from utils_entrezpy import MyMultiThreading_ncbi



import logging
logg = logging.getLogger("nuccoredb_retrival")


##################################################
# FUNCTIONS
##################################################

class MyMultiThreading_ncbi_mod(MyMultiThreading_ncbi):
  """
  Add timer and debugger to `run_batches_files` function 
  """

  @timer(my_logger=logg)
  def run_batches_files(self, ranges_batches):
    return super(MyMultiThreading_ncbi_mod, self).run_batches_files(ranges_batches)


##################################################
# ARGS
##################################################

def get_arg_parser():
    import argparse

    parser = argparse.ArgumentParser(description="Obtain IDs from NUCCORE DB",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--ncbi-api", required = True, 
                        help="NCBI API KEY")
    parser.add_argument("--query", required = True,
                        help="NCBI query" )
    parser.add_argument("-o", "--outfile", default = "sequences.fna", required = False,
                        help="Output file" )
    parser.add_argument("--log", default = "log.log", required = False)
    parser.add_argument("--threads", type = int, default = 1, required = False,
                        help="Threads to use" )
    
    return parser


##################################################
# MAIN
##################################################

if __name__ == "__main__":

    import os
    import glob
    from utils import run_cmd, split_by_size
    
    ARGS = get_arg_parser().parse_args()

    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "nuccoredb_retrival")
    my_logger = logger.get_logger()

    # Set NBCI API KEY
    os.environ["NCBI_API_KEY"] = f'{ARGS.ncbi_api}'

    # Obtain the number of results (nhits)
    cmd = f"esearch -db nuccore -query '{ARGS.query}' | xtract -pattern ENTREZ_DIRECT -element Count"
    my_logger.debug(f"Cmd to obtain nhits: {cmd}")
    cmd, nhits, p_stderr = run_cmd(cmd, split=False, shell=True)
    nhits = int(nhits)
    
    #
    ## Fetch the NUCCORE_ID in batches of 8000 elements 
    #

    ranges_batches = split_by_size(input=nhits, n=8000)
    ## Create batch queries
    my_logger.info(f"Downloading {nhits} NUCCORE UIDs in {len(ranges_batches)} batches using {ARGS.threads} threads...")
    file_tmp_prefix = "tmp.batch_0"
    args_query = [f"esearch -db nuccore -query '{ARGS.query}' | efetch -format docsum -start {batch[0]+1} -stop {batch[1]} | xtract -pattern DocumentSummary -element Id > {file_tmp_prefix}{i}" for i, batch in enumerate(ranges_batches)]
    args_cmd = [(query, False, True) for query in args_query]
    index_batches = {f"{file_tmp_prefix}{i}": i for i in range(len(ranges_batches))}
    ## Print the first and last 2
    [my_logger.debug(f"Args: {args_query[i]}") for i in range(2)]
    [my_logger.debug(f"Args: {args_query[i]}") for i in range(len(args_query)-2,len(args_query))]

    #
    ## Multiprocessing -> Get batches 
    #
    #NOTE: only 10 threads each iteration to prevent too many requests

    ncbi_batches = MyMultiThreading_ncbi_mod(args_cmd = args_cmd, myfunc = run_cmd, file_tmp = file_tmp_prefix, 
                          threads=ARGS.threads, my_logger = my_logger)
    ncbi_batches.run_batches_files(ranges_batches=ranges_batches)
    # Check correct downloading: in bytes and number of lines
    my_logger.info("Checking the size of the batch files...")
    last_batch_nlines= (ranges_batches[-1][1]) - ranges_batches[-1][0]
    ncbi_batches.check_batches_files(check_bytes= True, check_lines=True, 
                                    batch_nentries=8000, check_fasta=False,
                                    last_batch_nentries = last_batch_nlines, 
                                    index_batches=index_batches)
            
    # Join batches and eliminate them
    my_logger.info(f"Joinning data from batches...")
    ncbi_batches.join_batches_files(outfile = ARGS.outfile, header='NUCCORE_UID',
                                    sequentially=False)
    my_logger.info(f"File successfully generated: {ARGS.outfile}")

    
