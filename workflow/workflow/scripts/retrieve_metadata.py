#! /usr/bin/env/python
# coding: utf-8

## -----
# Retrival of NCBI data
# Author: G. Molano, LA (gonmola@hotmail.es)
# Last modified:
## -----

##################################################
# IMPORTS
##################################################

import csv # Do not remove
from utils_decorators import timer, debuger
from utils_my_logger import My_logger
from utils_entrezpy import MyMultiThreading_ncbi


import logging
logg = logging.getLogger("ncbi_request")

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

def parse_json_results(stdout, docsum_fun):
    """
    Parse json results from ncbi esummary function

    :params results: (list) STDOUT from esummary results
    :params docsum_fun: (func) Function to parse the json results
    :return (list) List of DocumentSummary objects
    """
    import json
    
    # Each data entry has to be splitted according to 'delim'
    delim = '{"header"'
    data = [delim+entry for res in stdout for entry in res.strip().split(delim) if entry]
    json_entries = [json.loads(entry.strip()) for entry in data]
    json_entries = [entry for entry in json_entries if 'result' in entry.keys()] # Filter error retrivals
    
    try:
      docsums = [docsum_fun(item) for entry in json_entries for key, item in entry['result'].items() if key != 'uids']
    except Exception as e:
        print(json_entries[0].keys())
        raise e
    
    return docsums

def elink_database_cmd(database, nuccore_uids):
    """
    Link nuccore UIDs with database uid
    :params database: (str) target database
    :params nuccore_uids: (list) list of nuccore UIDs (strings)

    :return (list) of target database uids
    """

    cmd_db_uids = f"""
    elink -db nuccore -id {','.join(nuccore_uids)} -target {database} -cmd neighbor | \
          xtract -pattern LinkSet -element Id | uniq
    """


    return cmd_db_uids

def esummary_cmd(database, database_uids):
    """
    Obtain docsums from ids
    :params database: (str) target database
    :params nuccore_uids: (list) list of nuccore UIDs (strings)

    :return (list) of target database uids
    """

    cmd_db_uids = f"esummary -id {','.join(database_uids)} -db {database} -mode json"


    return cmd_db_uids

def elink_database(database, nuccore_uids_all, ranges_batches, nhits, threads = 10, my_logger = None):
    from utils import run_cmd
    
    ## Multithreading -> Get batches #NOTE: only 10 threads each iteration to prevent too many requests

    # Get Linked UID
    my_logger.info(f"Linking {nhits} NUCCORE entries to {database} in {len(ranges_batches)} batches using {threads} threads...")
    args_query = [elink_database_cmd(database=database, nuccore_uids=nuccore_uids_all[start:end]) for start, end in ranges_batches]
    args_cmd = [(query, False, True) for query in args_query]

    ncbi_batches = MyMultiThreading_ncbi_mod(args_cmd = args_cmd, myfunc = run_cmd, 
                          threads=threads, my_logger = my_logger)
    results = ncbi_batches.run_batches_files(ranges_batches=ranges_batches)

    # Get data.frame between NUCCORE_UID and DB_UID
    data = [ row.split("\t") for res in results for row in res[1].strip().split('\n') if '\t']
    df_data = []
    for entry in data:
        if len(entry) > 2:
            entry_items = [[entry[0], item] for item in entry[1:]]
            df_data.extend(entry_items)
        if len(entry) == 2:
            df_data.append(entry)
    df = pd.DataFrame(df_data, columns=['NUCCORE_UID', f'{database.upper()}_UID'])
    df = df.astype('int64')

    # Get Summary
    database_uids = [str(item) for item in set(df[f'{database.upper()}_UID']) if isinstance(item, int)]
    ranges_batches = split_by_size(input=len(database_uids), n=ARGS.batch_size)
    my_logger.info(f"Fetching {len(database_uids)} {database} entries in {len(ranges_batches)} batches using {threads} threads...")
    
    args_query = [esummary_cmd(database=database, database_uids=database_uids[start:end]) for start, end in ranges_batches]
    args_cmd = [(query, False, True) for query in args_query]

    ncbi_batches = MyMultiThreading_ncbi_mod(args_cmd = args_cmd, myfunc = run_cmd, 
                          threads=threads, my_logger = my_logger)
    results = ncbi_batches.run_batches_files(ranges_batches=ranges_batches)


    # Get summary from UID

    return df, results
##################################################
# ARGS
##################################################
def get_args():
    import argparse

    parser = argparse.ArgumentParser(description="Obtain plasmid NCBI query ",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Mandatory argument
    parser.add_argument("--esearch-query", required = True, 
                        help="Terms to include in the esearch query")
    parser.add_argument("--ncbi-api", required = True, 
                        help="NCBI API KEY")
    parser.add_argument("--email", required = True,
                        help="Valid email for NCBI query" )
    parser.add_argument("--batch-size", type = int, default = 8000, required = False,
                        help="Number of elements in each batch" )
    parser.add_argument("-o", "--outdir", default = "results", required = False,
                        help="Output directory" )
    parser.add_argument("--log", default = "log.log", required = False)
    parser.add_argument("--threads", type = int, default = 1, required = False,
                        help="Threads to use" )
    return parser

##################################################
# MAIN
##################################################

if __name__ == '__main__':
    
    import re
    import pandas as pd
    import os
    from os.path import join
    from utils import run_cmd, split_by_size
    from utils_entrezpy import Docsum_nuccore, Docsum_assembly, Docsum_biosample
    import json

    ARGS = get_args().parse_args()
  
    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "ncbi_request")
    my_logger = logger.get_logger()

    # Set NBCI API KEY
    os.environ["NCBI_API_KEY"] = f'{ARGS.ncbi_api}'

    #
    ## NUCCORE
    #

    # Obtain the number of results (nhits)
    cmd = f"esearch -db nuccore -query '{ARGS.esearch_query}' | xtract -pattern ENTREZ_DIRECT -element Count"
    my_logger.debug(f"Cmd to obtain nhits: {cmd}")
    cmd, nhits, p_stderr = run_cmd(cmd, split=False, shell=True)
    nhits = int(nhits)
    
    #
    ## Fetch the HITS in batches of ARGS.batch_size elements 
    #

    ranges_batches = split_by_size(input=nhits, n=ARGS.batch_size)
    ## Create batch queries
    my_logger.info(f"Downloading {nhits} NUCCORE entries in {len(ranges_batches)} batches using {ARGS.threads} threads...")
    args_query = [f"esearch -db nuccore -query '{ARGS.esearch_query}' | esummary -mode json -start {start+1} -stop {end}" for start, end in ranges_batches]
    args_cmd = [(query, False, True) for query in args_query]
    ## Print the first and last 2
    [my_logger.debug(f"Args: {args_query[i]}") for i in range(2)]
    [my_logger.debug(f"Args: {args_query[i]}") for i in range(len(args_query)-2,len(args_query))]

    #NOTE: only 10 threads each iteration to prevent too many requests
    ncbi_batches = MyMultiThreading_ncbi_mod(
       args_cmd = args_cmd, myfunc = run_cmd,
       threads=ARGS.threads, my_logger = my_logger)
    results = ncbi_batches.run_batches_files(ranges_batches=ranges_batches)

    # Parse results in json format.
    stdout_res = [res[1] for res in results]
    docsums = parse_json_results(stdout=stdout_res, docsum_fun=Docsum_nuccore)
    
    # Process results
    header_nuccore = ["NUCCORE_UID", "NUCCORE_ACC", "NUCCORE_Description", 
                "NUCCORE_CreateDate", "NUCCORE_Topology", "NUCCORE_Completeness", 
                "NUCCORE_TaxonID",  "NUCCORE_Genome", "NUCCORE_Length", 
                "NUCCORE_DuplicatedEntry","NUCCORE_Source", "NUCCORE_BiosampleID"]
    records_nuccore = [[i.uid, i.accessionversion, i.title,
                i.createdate, i.topology, i.completeness,
                i.taxid, i.genome, i.slen,
                i.assemblyacc, i.sourcedb, i.biosample] for i in docsums if i.error == False]
    

    # Create dataframe and download
    nuccore_df = pd.DataFrame(records_nuccore, columns = header_nuccore)
    nuccore_df.to_csv(join(ARGS.outdir, "nuccore_records.csv"), index=False)
    
    # Link the NUCCORE_UIDs in batches of ARGS.batch_size elements 
    nuccore_df = pd.read_csv(join(ARGS.outdir, "nuccore_records.csv"))
    nuccore_ids = [str(item) for item in set(nuccore_df['NUCCORE_UID']) if isinstance(item, int)]
    nhits = len(nuccore_ids)
    ranges_batches = split_by_size(input=nhits, n=ARGS.batch_size)
    
    #
    ## ASSEMBLY
    # 
    
    # Create batch queries and execute
    assembly_linked_df, results = elink_database(database='assembly', nuccore_uids_all = nuccore_ids, 
                   ranges_batches = ranges_batches, nhits = nhits, 
                   threads = ARGS.threads, my_logger = my_logger)
    # Downloand table NUCCORE_UID - DB_UID
    assembly_linked_df.to_csv(join(ARGS.outdir, "assembly_linked_records.csv"), index=False)

    # Parse results in json format.
    stdout_res = [res[1] for res in results]
    docsums = parse_json_results(stdout=stdout_res, docsum_fun=Docsum_assembly)

    header_assembly = [
      "ASSEMBLY_UID", "ASSEMBLY_ACC", "ASSEMBLY_Status",
      "ASSEMBLY_coverage",
      "ASSEMBLY_SeqReleaseDate", "ASSEMBLY_SubmissionDate", 
      "ASSEMBLY_Lastest", "ASSEMBLY_BiosampleID"]
    # Keep only lastest assembly records
    records_assembly = [[
      i.uid, i.assemblyaccession, i.assemblystatus,
      i.coverage, 
      i.seqreleasedate,i.submissiondate, 
      True, i.biosample 
      ] for i in docsums if i.error == False if 'latest' in i.propertylist]

    # Create dataframe and download
    assembly_df = pd.DataFrame(records_assembly, columns = header_assembly)
    assembly_df = assembly_df.drop_duplicates()
    my_logger.info(f"Linked {len(assembly_df.index)} assembly entries")
    assembly_df.to_csv(join(ARGS.outdir, "assembly_records.csv"), index=False)


    #
    ## Biosample
    # 
    
    # Create batch queries and execute
    biosample_linked_df, results = elink_database(database='biosample', nuccore_uids_all = nuccore_ids, 
                   ranges_batches = ranges_batches, nhits = nhits, 
                   threads = ARGS.threads, my_logger = my_logger)
    # Downloand table NUCCORE_UID - DB_UID
    biosample_linked_df.to_csv(join(ARGS.outdir, "biosample_linked_records.csv"), index=False)

    # Parse results in json format.
    stdout_res = [res[1] for res in results]
    docsums = parse_json_results(stdout=stdout_res, docsum_fun=Docsum_biosample)

    header_biosample = [
      "BIOSAMPLE_UID", "BIOSAMPLE_ACC", "BIOSAMPLE_Location", 
      "BIOSAMPLE_Coordinates", "BIOSAMPLE_IsolationSource", "BIOSAMPLE_Host", 
      "BIOSAMPLE_CollectionDate","BIOSAMPLE_HostDisease", "BIOSAMPLE_SampleType"]
    records_biosample = [[
      i.uid, i.accession, i.sampledata.geographiclocation,
      i.sampledata.coordinates, i.sampledata.isolationsource, i.sampledata.host,
      i.sampledata.collectiondate, i.sampledata.hostdisease, i.sampledata.sampletype
      ] for i in docsums if i.error == False]
    
    # Create dataframe and download
    biosample_df = pd.DataFrame(records_biosample, columns = header_biosample)
    biosample_df = biosample_df.drop_duplicates()
    my_logger.info(f"Linked {len(biosample_df.index)} biosample entries")
    biosample_df.to_csv(join(ARGS.outdir, "biosample_records.csv"), index=False)


    # Join Nuccore, Assembly and Biosample records by "NUCCORE_UID"

    # Keep only nuccore_df entries bq are the target ones
    records_df = pd.merge(nuccore_df, biosample_df, how="left", 
                            left_on="NUCCORE_BiosampleID", right_on="BIOSAMPLE_ACC")
    # Keep only the previous merged entries bq are the target ones
    records_df = pd.merge(records_df, assembly_linked_df.reset_index(drop=True), how = "left",
                          left_on="NUCCORE_UID", right_on="NUCCORE_UID")
    # Keep only the assembly_df entries bq are the 'lastest' assemblies 
    records_df = pd.merge(records_df, assembly_df, how = "right",
                      left_on=["ASSEMBLY_UID"], right_on=["ASSEMBLY_UID"])
    
    records_df.to_csv(join(ARGS.outdir, "NAB_records.csv"), index=False)

    
  

