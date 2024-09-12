##################################################
# IMPORTS
##################################################
from utilsmeta.my_logger import My_logger
from utilsmeta import utils, parse_ete4
import utilsmeta
import pandas as pd
import numpy as np
from tqdm import tqdm
from copy import copy

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    
    # Logging
    logger = My_logger(log_filename = snakemake.log[0], logger_name = "infer_host")
    my_logger = logger.get_logger()
    
    ETE4 = parse_ete4.NCBITaxa_mod()

    #
    ## Here starts
    # 

    bio = utils.load_table(snakemake.input.bio,table="BIOSAMPLE")
    # BIOSAMPLE vars with sample source information
    host_vars = utilsmeta.ncbi.biosample_host_vars()
    bio[host_vars['taxid'][0]] = [i if isinstance(i, int) else -1 for i in list(bio[host_vars['taxid'][0]].fillna(-1)) ]

    # Previous curated mapping
    df = pd.read_csv(str(snakemake.params.mapping_checked))
    df['TAXONOMY_UID'] = [i for i in df['TAXONOMY_UID'].fillna("-1")]
    inmapping = utils.create_mapping(df,
        input_cols=['query'],
        output_cols = ['tags', 'TAXONOMY_UID'])

    # Load plasmid + Unify missing information labels
    names = []
    for i,row in bio.iterrows():
        names_row = [row[var] for var in host_vars['name'] if row[var]]
        names_row.extend([row[var] for var in host_vars['type'] if row[var]])
        if row[host_vars['taxid'][0]] != -1:
            names_row.append(f"taxid:{row[host_vars['taxid'][0]]}")
        
        names.append([i for i in names_row if isinstance(i,str)])

    bio['ECOSYSTEM_query'] = [','.join(i) for i in names]

    eco_map = []
    not_found = set()
    for i in tqdm(names):
        # For each record, join terms
        idict = {'tags' : set(), 'TAXONOMY_UID': set()}
        if i:
            for j in i:
                if 'taxid' in j:
                    match = {'tags':'host_associated','TAXONOMY_UID': int(j.replace('taxid:', ''))}
                else:
                    try:
                        match = inmapping[j]
                    except KeyError:
                        try:
                            match = inmapping[utils.process_string(j)]
                        except KeyError as e:
                            not_found.add(j)
                                
                # Aggregate terms
                idict['tags'].update(str(match['tags']).split(','))
                if match['TAXONOMY_UID'] != "-1":
                    try:
                        taxid = int(match['TAXONOMY_UID'])
                        idict['TAXONOMY_UID'].add(taxid)
                    except ValueError:
                        l = [int(i) for i in match['TAXONOMY_UID'].split(',')]
                        for i in l:
                            idict['TAXONOMY_UID'].add(i)
            
            idict['tags'] = [i for i in sorted(list(idict['tags'])) if i if i!='nan']

            # If various taxids, check wether is already contained in the lineage
            if len(idict['TAXONOMY_UID']) > 1:
                # Retrieve lineage of each taxid
                IDs = list(idict['TAXONOMY_UID'])
                lineage_ids = [ETE4.get_lineage(taxid) for taxid in IDs]

                # Check if taxid contained in other lineage
                eval = []
                for i, taxid in enumerate(IDs):
                    lineage_ids_copy = copy(lineage_ids)
                    lineage_ids_copy.pop(i)
                    bools = sum([True if taxid in lineage else False for lineage in lineage_ids_copy ])
                    eval.append(bools)

                # If at least one contained in other lineage
                if max(eval) > 0:
                    # Retrieve the not contained (bools == 0)
                    index_list = [i for i, item in enumerate(eval) if item == min(eval)]
                    idict['TAXONOMY_UID'] = [IDs[index_list[0]]]
                    
            eco_map.append(idict)
        else:
            eco_map.append(idict)
    
    if not_found:
        print(len(not_found))
        print(sorted(list(not_found)))
        # print('\n'.join(sorted(list(not_found))))
        raise KeyError
    
    bio['ECOSYSTEM_tags'] =  [','.join(i['tags']) for i in eco_map]
    bio['ECOSYSTEM_taxid'] =  [list(i['TAXONOMY_UID']) for i in eco_map]
    
    taxid_str = []
    for i in list(bio['ECOSYSTEM_taxid']):
        if len(i) == 0:
            ID = -1
        elif len(i) == 1:
            ID = i[0]
        else:
            ID = ','.join([str(j) for j in i])
        taxid_str.append(ID)
    bio['ECOSYSTEM_taxid'] = taxid_str

    # Flag non-uniques TAXONOMY_UIDS
    non_uniques = [True if len(item['TAXONOMY_UID']) > 1 else False for item in eco_map]
    bio['ECOSYSTEM_check'] = False
    bio.loc[non_uniques, 'ECOSYSTEM_check'] = True

    # Check previous infer mapping
    previous_mapping = pd.read_csv(snakemake.params.previous_mapping).set_index('BIOSAMPLE_UID')
    
    for index, row in bio.loc[bio['ECOSYSTEM_check']==True].iterrows():
        uid = row['BIOSAMPLE_UID']
        if row['ECOSYSTEM_query'] in list(previous_mapping['ECOSYSTEM_query']):
            bio.loc[bio['BIOSAMPLE_UID']==uid, "ECOSYSTEM_tags"] = previous_mapping.loc[uid,"ECOSYSTEM_tags"]
            bio.loc[bio['BIOSAMPLE_UID']==uid, "ECOSYSTEM_taxid"] = previous_mapping.loc[uid,"ECOSYSTEM_taxid"]
            bio.loc[bio['BIOSAMPLE_UID']==uid, "ECOSYSTEM_check"] = False
    
    # Save
    bio = bio.sort_values(by=['ECOSYSTEM_query', 'ECOSYSTEM_tags'])
    bio = bio.filter(regex='ECOSYSTEM|BIOSAMPLE_UID')

    bio.to_csv(snakemake.output[0], index=False)
    
    