import argparse
from utils_my_logger import My_logger
import os
from utils import run_cmd


##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--databases', '-db', help="The databases to be updated.", required=True, nargs='+')
    parser.add_argument('--outfile', '-o', help="Output file.", required=True)
    parser.add_argument("--log", default = "log.log", required = False)
    return parser


##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    from os.path import join
    from os import environ

    ARGS = get_arg_parser().parse_args()
    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "ncbi_request")
    my_logger = logger.get_logger()

    for db in ARGS.databases:
        # run ABRicate to update the database
        cmd = f"abricate-get_db --db {db} --force"
        my_logger.info(f'Update with ABRicate DB {db}: {cmd}')
        cmd, stdout, stderr = run_cmd(cmd, split=True, shell=False)
        
        with open(ARGS.outfile, 'a') as outfile:
            outfile.write(stdout)

        # reformat FASTA: remove non-UTF-8 characters
        input = join(environ["CONDA_PREFIX"], 'db', db, 'sequences')
        cmd = f"mv {input} {input}.tmp && iconv -c -t UTF-8 < {input}.tmp > {input}"
        my_logger.info(f'Remove non-UTF-8 symbols: {cmd}')
        run_cmd(cmd, split=False, shell=True)

    with open(ARGS.outfile, 'a') as outfile:
        outfile.write('\nUPDATE DONE.\n')