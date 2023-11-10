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
import numpy as np
from utils import create_mapping, unify_missing_info
from process_create_host_mapping import fix_host

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
    parser.add_argument("--input-pls", help="csv file for plasmid metadata")
    parser.add_argument("--mapping", help="csv file for finding BIOSAMPLE_Host")
    parser.add_argument("--regexes", help="Regex to find environmental samples",
                        type=json.loads )
    parser.add_argument("--log")
    parser.add_argument("--outfile")

    return parser

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    ARGS = get_parser().parse_args()
    
    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "infer_host")
    my_logger = logger.get_logger()

    #
    ## Here starts
    # 

    # Load mappings
    mapping_df = pd.read_csv(ARGS.mapping)
    map  = create_mapping(mapping_df, input_column="Input_host",
                   output_column="Final_host")

    # Load plasmid + Unify missing information labels
    df_pls = pd.read_csv(ARGS.input_pls)
    df_pls["BIOSAMPLE_Host"] = df_pls["BIOSAMPLE_Host"].astype("string")
    df_pls["BIOSAMPLE_IsolationSource"] = df_pls["BIOSAMPLE_IsolationSource"].astype("string")
    my_logger.info(f"Plasmid df: {len(df_pls.index)} entries")
    
    df_pls = unify_missing_info(df_pls, value_missing="nan", my_logger=my_logger)
    host_names = list(df_pls["BIOSAMPLE_Host"])
    iso_names = list(df_pls["BIOSAMPLE_IsolationSource"])
    
    # Assign hostname using host
    my_logger.info(f"Assigning host using BIOSAMPLE_Host information...")
    bio_assign_results = [(fix_host(
        host = name, mapping = map, regexes = ARGS.regexes,
        ask_ncbi = False, ncbi_db = None, 
        ask_two_words = True), "BIOSAMPLE_Host") for name in host_names]
    
    assigned_bio = [(res[1], source) for res, source in bio_assign_results if (res[1] != "unknown") and (res[1] != False) ]
    my_logger.info(f"Host_bio_assigned = {len(assigned_bio)}")
    my_logger.info(f"Unassigned = {len(bio_assign_results) - len(assigned_bio)}")
    
    # Assign host name using IsolationSource
    my_logger.info(f"Assigning remaining hosts using BIOSAMPLE_IsolationSource information...")
    assign_results = [(
        fix_host(
            host = iso, mapping = map, regexes = ARGS.regexes,
            ask_ncbi = False, ncbi_db = None, ask_two_words = True) ,
        "BIOSAMPLE_IsolationSource"
        ) if (bio[0][1] == "unknown") or (bio[0][1] == False) else bio for bio, iso in zip(
            bio_assign_results, iso_names)]
    
    assigned_iso = [(res, source)for res, source in assign_results if (res[1] != "unknown") and (res[1] != False) and(source == "BIOSAMPLE_IsolationSource")]
    my_logger.info(f"Host_iso_assigned = {len(assigned_iso)}")
    assigned_host = [(res, source) for res, source in assign_results if (res[1] != "unknown") and (res[1] != False)]
    my_logger.info(f"Total Host_assigned = {len(assigned_host)}")
    my_logger.info(f"NOT Host_assigned = {len(host_names) - len(assigned_host)}")
    
    # Add Host processed to plasmid information [Converting False to unknown]
    host_processed = [(res[1], source) if res[1] != False else ("unknown", source) for res, source in assign_results]
    df = pd.DataFrame(host_processed, columns=["BIOSAMPLE_Host_processed", "BIOSAMPLE_Host_processed_source"])
    df_pls = pd.concat([df_pls, df], axis=1)

    # Get final Host label
    # If
    special_category = ["unknown", "environmental"]
    df_pls["BIOSAMPLE_Host_label"] = np.where(
        # If Label is unknown or environmental
        df_pls["BIOSAMPLE_Host_processed"].isin(special_category),
            np.where(
                # If Assigned source is BIOSAMPLE_Host
                df_pls["BIOSAMPLE_Host_processed_source"] == "BIOSAMPLE_Host",
                    # Use the assigned name + (original name)
                    df_pls["BIOSAMPLE_Host_processed"] + " (" + df_pls["BIOSAMPLE_Host"] + ")",
                    # If the source is IsolationSource, check if the column is empty
                    np.where(
                        df_pls["BIOSAMPLE_IsolationSource"] == "nan",
                             df_pls["BIOSAMPLE_Host_processed"] + " (" + df_pls["BIOSAMPLE_Host"] + ")",
                             df_pls["BIOSAMPLE_Host_processed"] + " (" + df_pls["BIOSAMPLE_IsolationSource"] + ")"
                             )
                    ),
            # Use the assigned name
            df_pls["BIOSAMPLE_Host_processed"]
    )

    # Save
    df_pls.to_csv(ARGS.outfile, index=False)
    
    