#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 
# Authors:   Alejandra Gonzalez (gonmola@hotmail.es)
# ---

##################################################
# IMPORTS
##################################################


from utils_my_logger import My_logger
import pandas as pd
import os
from utils import run_cmd
from utils_entrezpy import MyMultiThreading_ncbi
from utils_decorators import timer

import logging
logg = logging.getLogger("fetch_fasta")

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
    
    parser = argparse.ArgumentParser(description="Fecth Fasta sequences from NUCCORE DB",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input-file", required = True, 
                        help="CSV file containning NUCCORE_UID column")
    parser.add_argument("--ncbi-api", required = True, 
                        help="NCBI API KEY")
    parser.add_argument("--email", required = True,
                        help="Valid email for NCBI query" )
    parser.add_argument("--batch-size", type = int, default = 8000, required = False,
                        help="Number of elements in each batch" )
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

    from utils import split_by_size

    ARGS = get_arg_parser().parse_args()
    
    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "fetch_fasta")
    my_logger = logger.get_logger()
    
    # Set NBCI API KEY
    os.environ["NCBI_API_KEY"] = f'{ARGS.ncbi_api}'

    # Parse input file and obtain list NUCCORE_UID
    df = pd.read_csv(ARGS.input_file)
    nuccore_ids = list(df['NUCCORE_UID'])
    
    #
    ## Fetch the NUCCORE_ID in batches of ARGS.batch_size elements 
    #
    ranges_batches = split_by_size(input=len(nuccore_ids), n=ARGS.batch_size)
    nuccore_ids_batch = [nuccore_ids[start:end] for start, end in ranges_batches]
    ## Create batch queries
    my_logger.info(f"Downloading {len(nuccore_ids)} NUCCORE UIDs in {len(ranges_batches)} batches using {ARGS.threads} threads...")
    file_tmp_prefix = "tmp.batch_0"
    args_query = [f"efetch -db nuccore -id {','.join([str(ids) for ids in batch])} -format fasta >  {file_tmp_prefix}{i}" for i, batch in enumerate(nuccore_ids_batch)]
    args_cmd = [(query, False, True) for query in args_query]
    index_batches = {f"{file_tmp_prefix}{i}": i for i in range(len(nuccore_ids_batch))}
    # print(index_batches)
    
    #
    ## Multithreading -> Get batches 
    #
    #NOTE: only 10 threads each iteration to prevent too many requests
    ncbi_batches = MyMultiThreading_ncbi_mod(args_cmd = args_cmd, myfunc = run_cmd, 
                                             file_tmp = file_tmp_prefix, threads=ARGS.threads, 
                                             my_logger = my_logger)
    ncbi_batches.run_batches_files(ranges_batches=ranges_batches)
    # Check correct downloading: in bytes and number of lines
    my_logger.info("Checking the size of the batch files...")
    last_batch_nentries= (ranges_batches[-1][1]) - ranges_batches[-1][0]
    ncbi_batches.check_batches_files(check_bytes= True, check_lines=False,
                                     check_fasta=True,
                                     batch_nentries=ARGS.batch_size, 
                                     last_batch_nentries=last_batch_nentries,
                                     index_batches=index_batches)
            
    # Join batches and eliminate them
    my_logger.info(f"Joinning data from batches...")
    ncbi_batches.join_batches_files(outfile = ARGS.outfile, header=None, 
                                    sequentially=True, compress=False)
    my_logger.info(f"File successfully generated: {ARGS.outfile}")

    
