##################################################
# IMPORTS
##################################################
from utilsmeta.my_logger import My_logger
from utilsmeta import utils
import pandas as pd
import numpy as np
from tqdm import tqdm


# Logging
logger = My_logger(log_filename = snakemake.log[0], logger_name = "infer_host")
my_logger = logger.get_logger()

#
## Here starts
# 

bio_loc = utils.load_table(snakemake.input.bio_loc,table="BIOSAMPLE").reset_index(drop=True)
eco_dis = pd.read_csv(snakemake.params.eco_dis_infer)
eco_dis['BIOSAMPLE_UID'] = eco_dis['BIOSAMPLE_UID'].astype(int) 

bio = bio_loc.merge(eco_dis, on="BIOSAMPLE_UID", how='left')
bio.set_index('BIOSAMPLE_UID', drop=False, inplace=True, verify_integrity=True)

# Load correction for locations
loc_corrections = pd.read_csv(snakemake.params.loc_corrections, dtype=str, na_values=[''])
correction_dict = utils.create_mapping(loc_corrections, 
                input_cols=["query_input"], 
                output_cols=["query_output", "notes"])


# Flag non-disease info with ECOSYSTEM_tags containing disease and viceversa: TO DO
bio['LOCATION_ecotag'] = [True if pd.notna(i) and 'location' in i.split(',') else False for i in bio['ECOSYSTEM_tags']]
bio['LOCATION_diseasetag'] = [True if pd.notna(i) and 'location' in i.split(',') else False for i in bio['DISEASE_tags']]

bio['LOCATION_check'] = np.where((
    (bio['LOCATION_ecotag']==True) | (bio['LOCATION_diseasetag']==True)
    ) & (pd.isna(bio['LOCATION_lat'])), True, False)

bio = bio.sort_values(by=['LOCATION_query', 'ECOSYSTEM_query', 'DISEASE_query'])

if bio['LOCATION_check'].any():
    df = bio[bio['LOCATION_check']==True]
    my_logger.info(f"To check: {len(df.index)}")

    # LOCATION queries
    queries = np.where(df['LOCATION_ecotag'] == True, df['ECOSYSTEM_query'], df['DISEASE_query'])

    # Add corresponding LOCATION query
    bio.loc[bio['LOCATION_check']==True, "LOCATION_query"] = np.where(df['LOCATION_ecotag'] == True, df['ECOSYSTEM_query'], df['DISEASE_query'])
    bio = bio.filter(regex='BIOSAMPLE_UID|LOCATION|DISEASE|ECOSYSTEM')

else:
    my_logger.info("Nothing to check: No location in ECOSYSTEM or DISEASE tags without location info")
    bio = bio.filter(regex='BIOSAMPLE_UID|LOCATION')

bio.to_csv(snakemake.output[0], index=False)


    

