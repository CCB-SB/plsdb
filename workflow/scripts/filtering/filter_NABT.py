#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 
# Authors:   Alejandra Gonzalez (gonmola@hotmail.es)
# ---

##################################################
# IMPORTS
##################################################

import pandas as pd
from utilsmeta.my_logger import My_logger
from utilsmeta.utils import load_table


##################################################
# MAIN
##################################################

if __name__ == "__main__":

    # Logging
    logger = My_logger(log_filename = snakemake.log[0], 
                       logger_name = "NABT_filtering")
    my_logger = logger.get_logger()

    # READ
    nucc = load_table(snakemake.input.nucc, table="NUCCORE")
    df_filter = pd.read_csv(snakemake.input.df_filt, low_memory=False)
                
    
    # Filter metadata and Save
    md_dict = {'nucc': 'NUCCORE_UID',
               'bio': 'BIOSAMPLE_UID',
               'ass': 'ASSEMBLY_UID',
               'tax': 'TAXONOMY_taxon_lineage'}

    for key, df_col in md_dict.items():
        if key in dfs:
            my_logger.info(f"Metadata {key} - {len(dfs[key].index)} entries")

            if key == 'nucc':
                dfs[key] = dfs[key][dfs[key][df_col].isin(list(df_filter[df_col]))]
            else:
                dfs[key] = dfs[key][dfs[key][df_col].isin(set(df_filter[df_col]))]
            
            my_logger.info(f"Filtered {key} - {len(dfs[key].index)} entries")
            
            dfs[key].to_csv(snakemake.output[key], index=False)