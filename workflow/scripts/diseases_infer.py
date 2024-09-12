##################################################
# IMPORTS
##################################################
from utilsmeta.my_logger import My_logger
from utilsmeta import utils
from rapidfuzz import fuzz
from utilsmeta.parse_ontology import Disease_ontology, Ontology
import utilsmeta
import pandas as pd
import numpy as np
from tqdm import tqdm


# Logging
logger = My_logger(log_filename = snakemake.log[0], logger_name = "infer_host")
my_logger = logger.get_logger()

# Load Disease Ontology information
DO = Disease_ontology(location=snakemake.params.disease_ont, name="DiseaseOntology")
SYMP = Ontology(location=snakemake.params.symp_ont, name="SymptomOntology")
#
## Here starts
# 

bio = utils.load_table(snakemake.input.bio,table="BIOSAMPLE")
eco = pd.read_csv(snakemake.params.eco_mapping)
eco['BIOSAMPLE_UID'] = eco['BIOSAMPLE_UID'].astype(int) 
bio = bio.reset_index(drop=True).merge(eco, on="BIOSAMPLE_UID", how='left')

# BIOSAMPLE vars with sample source information
disease_vars = utilsmeta.ncbi.biosample_host_vars()['health_status']

# Previous curated mapping
df = pd.read_csv(str(snakemake.params.mapping_checked))
df['ONTOLOGY_ID'] = [i for i in df['ONTOLOGY_ID'].fillna(-1)]
inmapping = utils.create_mapping(df,
    input_cols=['query'],
    output_cols = ['tags', 'ONTOLOGY_ID'])

# Load plasmid + Unify missing information labels
names = []
for i,row in bio.iterrows():
    names_row = [row[var] for var in disease_vars if row[var]]
    names.append([i for i in names_row if isinstance(i,str) if not i.isnumeric()])

bio['DISEASE_query'] = [','.join(i) for i in names]

eco_map = []
not_found = set()
for i in tqdm(names):
    # For each record, join terms
    idict = {'tags' : set(), 'ONTOLOGY_ID': set()}
    if i:
        for j in i:
            try:
                match = inmapping[j]
            except KeyError:
                try:
                    match = inmapping[utils.process_string(j)]
                except KeyError as e:
                    not_found.add(utils.process_string(j))
                            
            # Aggregate terms
            idict['tags'].update(str(match['tags']).split(','))
            if match['ONTOLOGY_ID'] != "-1":
                l = [i for i in match['ONTOLOGY_ID'].split(',')]
                for i in l:
                    idict['ONTOLOGY_ID'].add(i)
    
        idict['tags'] = [i for i in sorted(list(idict['tags'])) if i if i!='nan']
                
        eco_map.append(idict)
    else:
        eco_map.append(idict)

if not_found:
    print(len(not_found))
    print(not_found)
    # print(sorted(list(not_found)))
    # print('\n'.join(sorted(list(not_found))))
    raise KeyError

bio['DISEASE_tags'] =  [','.join(i['tags']) for i in eco_map]
bio['DISEASE_ontid'] =  [list(i['ONTOLOGY_ID']) for i in eco_map]

# Ontid from list to str
ontid_str = []
for i in list(bio['DISEASE_ontid']):
    if len(i) == 0:
        ID = -1
    elif len(i) == 1:
        ID = i[0]
    else:
        ID = ','.join([str(j) for j in i])
    ontid_str.append(ID)
bio['DISEASE_ontid'] = ontid_str

# Flag non-disease info with ECOSYSTEM_tags containing disease and viceversa: TO DO
bio['DISEASE_ecotag'] = [True if pd.notna(i) and 'disease' in i.split(',') else False for i in bio['ECOSYSTEM_tags']]
bio['DISEASE_check'] = np.where((bio['DISEASE_ecotag']==True) & (bio['DISEASE_ontid']==-1), True, False)

bio['ECOSYSTEM_diseasetag'] = [True if pd.notna(i) and 'ecosystem' in i.split(',') else False for i in bio['DISEASE_tags']]
bio['ECOSYSTEM_check'] = np.where((bio['ECOSYSTEM_diseasetag']==True) & (bio['ECOSYSTEM_taxid']==-1), True, False)

bio = bio.sort_values(by=['DISEASE_query', 'DISEASE_tags', 'ECOSYSTEM_query', 'ECOSYSTEM_tags'])

# Check previous infer mapping
previous_mapping = pd.read_csv(snakemake.params.previous_mapping).set_index('BIOSAMPLE_UID')
    
for index, row in bio.loc[bio['DISEASE_check']==True].iterrows():
    uid = row['BIOSAMPLE_UID']
    if uid in list(previous_mapping.index):
        bio.loc[bio['BIOSAMPLE_UID']==uid, "DISEASE_tags"] = previous_mapping.loc[uid,"DISEASE_tags"]
        bio.loc[bio['BIOSAMPLE_UID']==uid, "DISEASE_ontid"] = previous_mapping.loc[uid,"DISEASE_ontid"]
        bio.loc[bio['BIOSAMPLE_UID']==uid, "DISEASE_check"] = False

if bio['DISEASE_check'].any():
    
    # Queries
    disease_names = [i for i in list(bio[bio['DISEASE_check']==True]["ECOSYSTEM_query"])]

    # Compute weighted ratio to DO
    wratio_results = [DO.get_matching_node(query=item, threshold = snakemake.params.threshold_wratio_do,
                                            scorer = fuzz.WRatio) for item in tqdm(disease_names)]
    matches_wratio = [i for i in wratio_results if i if i['ONTOLOGY_match'] != '']
    my_logger.info(f"Matches for weighted ratio (>={snakemake.params.threshold_wratio_do} similarity) = {len(matches_wratio)}")

    # Map to symptom ontology for unmatches results
    my_logger.info(f"Preprocessing and mapping disease names with <{snakemake.params.threshold_wratio_do} similarity in {DO.name} to {SYMP.name}")
    symp_results = [SYMP.get_matching_node(query=disease_names[i],
                                            threshold = snakemake.params.threshold_wratio_symp,
                                            scorer = fuzz.WRatio
                                            ) if not item or item['ONTOLOGY_match'] == '' or item['score'] < 95 else item for i, item in enumerate(wratio_results)]
    symp_results = [item if item and item['score'] >= 95 else wratio_results[i] for i, item in enumerate(symp_results)]
    matches_symp = [item for item in symp_results if item if item['ONTOLOGY_match'] != '' if item['ONTOLOGY_match'] == SYMP.name]
    my_logger.info(f"Total disease names match to a {SYMP.name} term = {len(matches_symp)} with (>= {snakemake.params.threshold_wratio_symp})")

    # Compute token set ratio for unmatches to DO
    tsr_results = [DO.get_matching_node(query=disease_names[i],
                                            threshold = snakemake.params.threshold_tsr_do,
                                            scorer=fuzz.token_set_ratio
                                                            ) if not item or item['ONTOLOGY_match'] == '' else item for i, item in tqdm(enumerate(symp_results))]
    matches = [i for i in tsr_results if i if i['ONTOLOGY_match'] != '']
    matches_tsr = [i for i in matches if i not in matches_wratio if i not in matches_symp]
    my_logger.info(f"Matches for token set ratio in {DO.name} (>={snakemake.params.threshold_tsr_do} similarity) = {len(matches_tsr)}")

    # Compute token set ratio for unmatches to SYMP
    final_results = [SYMP.get_matching_node(query=disease_names[i],
                                            threshold = snakemake.params.threshold_tsr_symp,
                                            scorer=fuzz.token_set_ratio
                                                            ) if not item or item['ONTOLOGY_match'] == '' else item for i, item in tqdm(enumerate(tsr_results))]
    matches = [i for i in final_results if i if i['ONTOLOGY_match'] != '']
    matches_tsr_symp = [i for i in matches if i not in matches_wratio if i not in matches_symp if i not in matches_tsr]
    my_logger.info(f"Matches for token set ratio in {SYMP.name} (>={snakemake.params.threshold_tsr_symp} similarity) = {len(matches_tsr_symp)}")

    my_logger.info(f"Total disease names match to a disease/symp ontology term = {len(matches)}")

    # Record mapping
    new_map = [] 
    for query, item in zip(disease_names, final_results):
        if item:
            new_map.append(item)
        else:
            new_map.append({'query': utilsmeta.utils.process_string(query), 'method': "not_assigned", "ONTOLOGY_ID": -1})
    new_mapping = pd.DataFrame.from_records(new_map)

    # Insert disease mapping
    new_mapping['DISEASE_query'] = "ECOSYSTEM_query"
    new_mapping['DISEASE_tags'] = new_mapping['tags']
    new_mapping['DISEASE_ontid'] = new_mapping['ONTOLOGY_ID']
    bio.loc[bio['DISEASE_check']==True, ["DISEASE_query", "DISEASE_tags", "DISEASE_ontid"]] = new_mapping[["DISEASE_query", "DISEASE_tags", "DISEASE_ontid"]].to_numpy().tolist()

# Append disease terms to enhance evaluation
terms = []
for i in tqdm(bio['DISEASE_ontid']):
    if i != -1 and i!="-1" and pd.notna(i):
        t = [DO.ont[j].name if 'DOID' in j else SYMP.ont[j].name for j in i.split(',')]
        t = ','.join(t)
    else:
        t=""
    terms.append(t)

# Insert disease group
disease_group = []
disease_group2 = []
disease_group3 = []
for i in tqdm(bio['DISEASE_ontid']):
    if i != -1 and i!="-1" and pd.notna(i):
        t_list = [list(DO.ont[j].superclasses()) if 'DOID' in j else list(SYMP.ont[j].superclasses()) for j in i.split(',')]
        t = ','.join([  j[-2].name for j in t_list if len(j) > 1])
        t2 = ','.join([ j[-3].name for j in t_list if len(j) > 2])
        t3 = ','.join([ j[-4].name for j in t_list if len(j) > 3])
    else:
        t=""
        t2=""
        t3=""
    disease_group.append(t)
    disease_group2.append(t2)
    disease_group3.append(t3)

bio = bio.filter(regex='DISEASE|ECOSYSTEM|BIOSAMPLE_UID')
bio.insert(8, "DISEASE_term", terms)
bio.insert(9, "DISEASE_group", disease_group)
bio.insert(10, "DISEASE_group2", disease_group2)
bio.insert(11, "DISEASE_group3", disease_group3)
bio.to_csv(snakemake.output[0], index=False)
    

