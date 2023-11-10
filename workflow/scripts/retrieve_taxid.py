#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 
# Authors:   Alejandra Gonzalez (gonmola@hotmail.es)
# ---

###########
# IMPORTS #
###########

from ete3 import NCBITaxa
import argparse
from utils_my_logger import My_logger


# NCBI db
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()


#############
# FUNCTIONS #
#############

def get_taxa_from_specie_taxid(sp_taxid, taxid=True, append_ID = False):
    """
    Takes a specie taxID and returns a list with all the taxonomy levels
    
    Arguments
    ----------
        sp_taxid: int or str
            Integer containing the specie taxID or str
        taxid: bool
            Wheter the sp_taxid is taxid (int) or str(description)
        append_ID: bool
            Wheter to append taxid to the final list
        
    
    Returns
    -------
    list
        Eight-lineage taxonomy (i.e., Kingdom, Phylum,..., Specie, Strain) 
    
    """
    import re
    # Define 8-level taxa + formating prefixes
    eight_lineage_list = ['superkingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species', 'strain']

    if not taxid:
        name2taxid = ncbi.get_name_translator([sp_taxid])
        if name2taxid:
            sp_taxid = list(name2taxid.values())[0][0]
        else:
            match = re.search(r'sp.?$', sp_taxid)
            if match:
                sp_taxid = sp_taxid.replace(match.group(), '').strip()
                return(get_taxa_from_specie_taxid(sp_taxid, taxid=False, 
                                           append_ID = True))
            else:
                raise ValueError(f"No taxid found for {sp_taxid}") 

    lineage = ncbi.get_lineage(sp_taxid) # contains the original order of taxids
    names = ncbi.get_taxid_translator(lineage)
    # Get only 8-lineage (taxids)
    ranks = ncbi.get_rank(lineage)
    # invert the ranks dictionary
    inv_ranks = {v: k for k, v in ranks.items()}

    value = []
    value_id = []
    last_name = ""
    last_rank = ""

    for rank in eight_lineage_list:
        if rank in inv_ranks:
            id = inv_ranks[rank]
            value.append(names[id])
            value_id.append(id)
            last_name = value
            last_rank = rank
        else:
            value.append("")
            value_id.append("")

    # replace spaces by _
    value = list(map(lambda x: x.replace(" ", "_"), value))

    if append_ID:
        value.append(sp_taxid)
        value_id.append(sp_taxid)

    return value, value_id, last_rank, last_name

if __name__ == "__main__":

    

    parser = argparse.ArgumentParser(description="Obtain taxid seven lineage",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Mandatory argument
    parser.add_argument("--input-file", required = True, 
                        help="CSV file containing merged records from nuccore, assembly, and biosample DBs.")
    # Optional args
    parser.add_argument("-o", "--outdir", default = "results", required = False,
                        help="Output directory" )
    parser.add_argument("--log", default = "log.log", required = False)


    args = parser.parse_args()
  
    # Logging
    logger = My_logger(log_filename = args.log, logger_name = "taxid_get_lineage")
    my_logger = logger.get_logger()

    #
    ## Here starts
    #
    import pandas as pd
    import numpy as np
    from os.path import join

    # Load file and obtain taxid list
    records_df = pd.read_csv(args.input_file)
    taxid_list = list(records_df["NUCCORE_TaxonID"])
    taxid_description = list(records_df["NUCCORE_Description"])

    # Get lineage
    eight_lineage_list = []
    eight_lineage_list_id = []
    taxon_rank = []
    last_name = []
    my_logger.info(f"Searching taxonomy lineage for {len(taxid_list)} taxids...")
    for i, ID in enumerate(taxid_list):
        try:
           lineage, lineage_id, last_rank, last_name_rank =  get_taxa_from_specie_taxid(ID, taxid=True, append_ID=True)
        except ValueError:
            # Search the two first fields of NUCCORE_Description
            description = ' '.join(taxid_description[i].split(" ")[:2])
            my_logger.info(f"NUCCORE_TaxonID ({ID}) not found...Extracting taxonomy information from NUCCORE_Description ({description})")
            lineage, lineage_id, last_rank, last_name_rank =  get_taxa_from_specie_taxid(description, taxid=False, 
                                                              append_ID=True)
        eight_lineage_list.append(lineage)
        taxon_rank.append(last_rank)
        last_name.append(last_name_rank)
        eight_lineage_list_id.append(lineage_id)
    
    # Taxonomy_names
    tax_df = pd.DataFrame(eight_lineage_list, 
                          columns=['TAXONOMY_superkingdom', 'TAXONOMY_phylum', 'TAXONOMY_class', 
                                   'TAXONOMY_order','TAXONOMY_family', 'TAXONOMY_genus', 
                                   'TAXONOMY_species', "TAXONOMY_strain", "TAXONOMY_UID"])
    tax_df['TAXONOMY_taxon_rank'] = taxon_rank
    tax_df['TAXONOMY_taxon_name'] = last_name
    tax_df["TAXONOMY_taxon_lineage"] =  tax_df[['TAXONOMY_superkingdom', 'TAXONOMY_phylum', 'TAXONOMY_class', 
                                   'TAXONOMY_order','TAXONOMY_family', 'TAXONOMY_genus', 
                                   'TAXONOMY_species', "TAXONOMY_strain"]].apply(lambda x: ';'.join(x.dropna()), axis=1)
    tax_df.to_csv(join(args.outdir, "taxonomy_records.csv"), index=False)
    # Taxonomy_Ids
    tax_df_id = pd.DataFrame(eight_lineage_list_id, 
                          columns=['TAXONOMY_superkingdom_id', 'TAXONOMY_phylum_id', 'TAXONOMY_class_id', 
                                   'TAXONOMY_order_id','TAXONOMY_family_id', 'TAXONOMY_genus_id', 
                                   'TAXONOMY_species_id', "TAXONOMY_strain_id", "TAXONOMY_UID"])
    tax_df_id.to_csv(join(args.outdir, "taxonomy_id_records.csv"), index=False)
    tax_df_id.drop(["TAXONOMY_UID"], axis = 1, inplace=True)
    
    # Merge DFs
    if 'TAXONOMY_superkingdom' in records_df.columns:
        records_df.drop(['TAXONOMY_superkingdom', 'TAXONOMY_phylum', 'TAXONOMY_class', 
                                    'TAXONOMY_order','TAXONOMY_family', 'TAXONOMY_genus', 
                                    'TAXONOMY_species', "TAXONOMY_strain", "TAXONOMY_UID"], axis=1, inplace=True)
    records_df = pd.concat([records_df, tax_df, tax_df_id], axis="columns")

    outfile = join(args.outdir, "NABT_records.csv")
    records_df.to_csv(outfile, index=False)
    
    my_logger.info(f"Successful retrival for {len(eight_lineage_list)} taxids...")
    my_logger.info(f"{outfile} created.")