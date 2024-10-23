from os.path import join, basename, splitext
import os
from glob import glob
from concurrent.futures import ThreadPoolExecutor
from utilsmeta.concurrency import split_by_size
from utilsmeta.my_logger import My_logger
from utilsmeta.utils import run_cmd

# Logging
logger = My_logger(log_filename = snakemake.log[0], logger_name = "ncbi_request")
my_logger = logger.get_logger()

os.makedirs(snakemake.output.DIR, exist_ok=True)
files = glob(join(snakemake.input.DIR, "*.gbk"))

arg_queries = [f"""ruby {snakemake.params.DIR_repo}cgview_builder_cli.rb --sequence {i} \
                --outfile {snakemake.output.DIR}/{splitext(basename(i))[0]}.json \
                --map_name name --config {snakemake.input.config}""" for i in files]
my_logger.info(f"CMD = {arg_queries[0]}")
arg_cmds = [(arg_query, False, True) for arg_query in arg_queries]

ranges_threads = split_by_size(input=len(arg_cmds), n=snakemake.threads)
my_logger.info(f"Ranges batches: len = {len(ranges_threads)}; [0] = {ranges_threads[0]}")

with ThreadPoolExecutor(snakemake.threads) as pool:
    try:
        # issue all tasks to the thread pool
        futures = [pool.submit(run_cmd, j[0], j[1]) for j in arg_cmds]
        # retrieve all return values in order
        results = [future.result() for future in futures]
    except Exception as e:
        my_logger.error(f"Something went wrong with Multithreading: {e}")

my_logger.info(f"Jobs done")
