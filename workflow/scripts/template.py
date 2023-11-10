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


##################################################
# FUNCTIONS
##################################################

##################################################
# ARGS
##################################################

def get_parser():
    import argparse
    import json
    # create parser object
    parser = argparse.ArgumentParser(
        description="Assign BIOSAMPLE Host using iso/host mappings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # specify our desired options 
    # by default ArgumentParser will add an help option 
    parser.add_argument("--input", help="csv file for plasmid metadata")
    parser.add_argument("--regexes", help="Regex to find environmental samples",
                        type=json.loads )
    parser.add_argument("--threads", type = int, default = 1, required = False,
        help="Threads to use" )
    parser.add_argument("--log")
    parser.add_argument("--outfile")

    return parser

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    ARGS = get_parser().parse_args()
    
    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "")
    my_logger = logger.get_logger()

    
    