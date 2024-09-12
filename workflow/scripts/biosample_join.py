selected_biovars = [
    'BIOSAMPLE_UID', 
    'NUCCORE_UID', 
    'BIOSAMPLE_ACC',
    'BIOSAMPLE_title', 
    'BIOSAMPLE_taxonomy', 
    'BIOSAMPLE_package', 
    'BIOSAMPLE_last_update',
    'BIOSAMPLE_metagenome_source',
    # Isolate
    'BIOSAMPLE_isolate', 
    'BIOSAMPLE_isolate_name_alias',
    'BIOSAMPLE_isol_growth_condt',
    'BIOSAMPLE_lab_host',
    'BIOSAMPLE_cell_type',
    # Plasmid Features
    'BIOSAMPLE_num_replicons',
    'BIOSAMPLE_extrachrom_elements',
    'BIOSAMPLE_encoded_traits',
    'BIOSAMPLE_propagation',
    # Plasmid Host features
    'BIOSAMPLE_genotype', 
    'BIOSAMPLE_biol_stat',
    'BIOSAMPLE_biovar', 
    'BIOSAMPLE_rel_to_oxygen', 
    'BIOSAMPLE_beta_lactamase_family',
    'BIOSAMPLE_carbapenemase',
    'BIOSAMPLE_edta_inhibitor_tested',
    'BIOSAMPLE_pathotype',
    'BIOSAMPLE_pathovar', 
    'BIOSAMPLE_pathogenicity',
    'BIOSAMPLE_phenotype',
    # Ecosystem HA relationship
    'BIOSAMPLE_host_dependence', 
    'BIOSAMPLE_biotic_regm',
    'BIOSAMPLE_biotic_relationship', 
    # Ecosystem HA Age 
    'BIOSAMPLE_host_age', 
    'BIOSAMPLE_age', 
    'BIOSAMPLE_host_life_stage',
    # Ecosystem HA Sex
    'BIOSAMPLE_host_sex', 
    'BIOSAMPLE_sex',
    # Ecosystem HA Disease state 
    'BIOSAMPLE_host_disease_stage', 
    'BIOSAMPLE_host_disease_outcome', 
    # Ecosystem HA Others
    'BIOSAMPLE_host_diet',
    'BIOSAMPLE_ethnicity',
    'BIOSAMPLE_plant_body_site',
    'BIOSAMPLE_cur_vegetation', 
    # Ecosystem Others
    'BIOSAMPLE_extreme_event',
    'BIOSAMPLE_extreme_salinity',
    ]

import pandas as pd
from utilsmeta import utils

bio = utils.load_table(snakemake.input.bio, table="BIOSAMPLE").reset_index(drop=True)
bio = bio.loc[:, selected_biovars]
dis_eco = pd.read_csv(snakemake.input.dis_eco)
loc = pd.read_csv(snakemake.input.loc)

bio = pd.merge(bio, dis_eco, how="left", on="BIOSAMPLE_UID")
bio = pd.merge(bio, loc, how="left", on="BIOSAMPLE_UID")

print(bio.head())
bio.to_csv(snakemake.output[0], index=False)

