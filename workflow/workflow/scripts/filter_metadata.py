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
from utils_decorators import timer, debuger

import logging
logg = logging.getLogger("plasmid_metadata_filtering")


##################################################
# FUNCTIONS
##################################################


@timer(my_logger=logg)
@debuger(my_logger=logg)
def filter_by_metadata(input_file, filter_regex, my_logger=None):
    """
        Filter plasmid based on metadata attributes:
            - Duplicated entry (NUCCORE_DuplicatedEntry)
            - Description (regex)
            - Assembly Lastest (True)
            - Assembly Completeness (Complete)
            - Taxonomy (Bacteria)
    """
    import pandas as pd
    import re

    def pls_filter(uid):

        assembly_tag = pls.loc[uid, 'ASSEMBLY_Status']

        has_assembly_tag = pd.notnull(pls.loc[uid, 'ASSEMBLY_UID']) and pd.notnull(assembly_tag) and assembly_tag != ""

        complete_tag  = pls.loc[uid, 'NUCCORE_Completeness']

        has_complete_tag = pd.notnull(complete_tag) and complete_tag != ""

        if has_assembly_tag and has_complete_tag:
            return assembly_tag == 'Complete Genome' and complete_tag == 'complete'
        elif not has_assembly_tag:
            return complete_tag == 'complete'
        elif not has_complete_tag:
            return assembly_tag == 'Complete Genome'
        return False


    # Read input file
    pls = pd.read_csv(input_file, header=0)
    pls = pls.drop_duplicates()
    pls.set_index('NUCCORE_UID', drop=False, inplace=True, verify_integrity=True)
    my_logger.info(f'Read in {len(pls.index)} plasmid records')

    # filter by Duplicated Entry
    nuccore_ids_dup = set(pls['NUCCORE_DuplicatedEntry'])
    pls = pls.loc[~pls['NUCCORE_ACC'].isin(nuccore_ids_dup),:]
    my_logger.info(f'Filtered: Duplicated entries (Refseq & insd, keep Refseq): kept {len(pls.index)} records')

    # filter by description
    keep_rows = pd.Series([re.search(filter_regex, d, flags=re.IGNORECASE) is None for d in pls['NUCCORE_Description']], index=pls.index)
    pls = pls.loc[keep_rows,:]
    my_logger.info(f'Filtered: description "{filter_regex}": kept {len(pls.index)} records')
    

    # filter by (assembly) completeness
    my_logger.info(f"Completeness nuccore tag: {set(pls['NUCCORE_Completeness'])}")
    my_logger.info(f"Completeness assembly tag: {set(pls['ASSEMBLY_Status'])}")
    pls = pls.loc[pls.index.map(pls_filter),:]
    my_logger.info(f'Filtered: (assembly) completeness: kept {len(pls.index)} records')

    # filter by taxonomy
    pls = pls.loc[pls['TAXONOMY_superkingdom'] == 'Bacteria',:]
    my_logger.info(f'Filtered: superkingdom ID is Bacteria: kept {len(pls.index)} records')

    return pls


##################################################
# ARGS
##################################################
def get_arg_parser():
    import argparse

    parser = argparse.ArgumentParser(description="Filter NCBI data",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Mandatory argument
    parser.add_argument("--input-file", required = True, 
                        help="CSVfile containing merged records from nuccore, assembly, biosample, and taxonomy DBs.")
    # Optional args
    parser.add_argument("--filter-regex", required=False, default=None,
                        help="Regex to filter NUCCORE_Description")
    parser.add_argument("-o", "--outfile", default = "md_filt.txt", required = False,
                        help="Output file" )
    parser.add_argument("--log", default = "log.log", required = False)

    return parser

##################################################
# MAIN
##################################################

if __name__ == "__main__":

    ARGS = get_arg_parser().parse_args()
    
    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "plasmid_metadata_filtering")
    my_logger = logger.get_logger()
    
    # Filter by metadata
    pls = filter_by_metadata(input_file = ARGS.input_file, filter_regex = ARGS.filter_regex,
                    my_logger = my_logger)

    pls.to_csv(ARGS.outfile, index=False)

