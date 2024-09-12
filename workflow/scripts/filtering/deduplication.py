##################################################
# FUNCTIONS
##################################################

def remove_duplicates(input_file, outfile, my_logger=None):
    import subprocess
    my_logger.info(f"Checking for possible duplicated sequences in {input_file}")

    tmp_file = "tmp.duplicated_ids.txt" # info of duplicated IDs
    cmd = f"cat {input_file} | seqkit rmdup --by-seq --dup-num-file {tmp_file} > {outfile}"
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    proc.wait()
    message = proc.stderr.read().decode()

    my_logger.info(f"Seqkit output: {message}")

    return tmp_file


##################################################
# MAIN
##################################################

if __name__ == "__main__":
    
    import numpy as np
    import pandas as pd
    from utilsmeta.utils import filter_fasta, merge_NABT, read_NABT
    from utilsmeta.my_logger import My_logger
    
    # Logging
    logger = My_logger(log_filename = snakemake.log[0], 
                       logger_name = "plasmid_sequence_filtering")
    my_logger = logger.get_logger()

    # Merge
    dfs = read_NABT(nucc_path=snakemake.input.nucc,
                    bio_path=snakemake.input.bio, 
                    ass_path=snakemake.input.ass,
                    tax_path=snakemake.input.tax)
    md = merge_NABT(nucc=dfs['nucc'], bio=dfs['bio'],
               ass=dfs['ass'], tax=dfs['tax'])
    print(md.head())

    ## Check if Assembly and Biosample Info
    criteria_dict = {
        'biosample': ['BIOSAMPLE_UID'],
        'assembly': ['ASSEMBLY_UID'],
        'location': [
            'BIOSAMPLE_geo_loc_name', 
            'BIOSAMPLE_lat_lon'],
        'host_status': [
            'BIOSAMPLE_host', 
            'BIOSAMPLE_host_disease',
            'BIOSAMPLE_host_health_state', 
            'BIOSAMPLE_phenotype'],
        'sample_source': [
            'BIOSAMPLE_isolation_source',
            'BIOSAMPLE_source_type'],
        'resistances': [
            'BIOSAMPLE_encoded_traits',
            'BIOSAMPLE_carbapenemase',
            'BIOSAMPLE_pathotype'],
        'host_extra': [
            'BIOSAMPLE_host_sex',
            'BIOSAMPLE_host_age'],
        'plasmid_extra': [
            'BIOSAMPLE_num_replicons']
    }
    for key, var_list in criteria_dict.items():
        for var in var_list:
            md[f'BOOL_{var}'] = np.where(md[var].isnull(), False, True)


    # Parse duplicated seqs id info in list of tuples
    dup_groups = [] 
    with open(snakemake.input.dup_ids, 'r') as f:
        for line in f:
            # [NUCCORE_ACC, NUCCORE_ACC,...,]
            nuccore_acc_list = line.strip().split("\t")[1].split(', ')
            dup_groups.append(nuccore_acc_list)
    ##
    # Define criteria
    ##

    selecting_criteria = """
    Compare two plasmid records, prefer:
    * RefSeq
    * Biosample entry
    * Biosample Location information
    * Assembly info
    * More recent Assembly release date
    * Biosample Host information
    * Biosample Host Disease information
    * Biosample Sample Source information
    * Biosample Resistances information
    * Biosample Host extra information (sex, age, ...,)
    * Biosample plasmid extra information (i.e., num_replicons)
    * Highest coverage
    * More recent Nuccore creation date
    * If all equals, choose the first one
    """
    sorting_columns = ['NUCCORE_Source']
    for i in ['biosample', 'location', 'assembly']:
        sorting_columns.extend([f"BOOL_{var}" for var in criteria_dict[i]])
    sorting_columns.append("ASSEMBLY_SeqReleaseDate")
    for i in ['host_status', 'sample_source', 'resistances', 'plasmid_extra']:
        sorting_columns.extend([f"BOOL_{var}" for var in criteria_dict[i]])
    sorting_columns.append("ASSEMBLY_coverage")
    sorting_columns.append("NUCCORE_CreateDate")

    my_logger.info(f"""There are {len(dup_groups)} groups of equal sequences.\nTo choose the nuccore record,
                   the criteria will be the following:\n {selecting_criteria}""")
    my_logger.info(f"Sorting columns: {sorting_columns}")

    # Compare records
    nequals = 0
    md_identical = pd.DataFrame()
    ids_to_delete = []
    for group in dup_groups:
        # subset md
        md_sub = md[md['NUCCORE_ACC'].isin(group)]
        
        # Sort by criteria
        md_sub = md_sub.sort_values(by=sorting_columns,
            ascending=False)

        # Check if first and second items are equals
        record_to_keep = md_sub.iloc[0]
        equals = record_to_keep[sorting_columns].equals(md_sub[sorting_columns].iloc[1])

        if equals:
            nequals = nequals + 1
        
        # Select the the first item and keep record of the others
        records_to_delete = md_sub[md_sub['NUCCORE_ACC']!= record_to_keep.loc['NUCCORE_ACC']]
        to_delete = list(records_to_delete.loc[:, ['NUCCORE_ACC']])
        ids_to_delete.extend(to_delete)

        # Keep metadata information from the duplicated item
        records_to_delete.loc[:,['NUCCORE_Identical']] = record_to_keep['NUCCORE_ACC']
        md_identical = pd.concat([md_identical, records_to_delete])
    
    my_logger.info(f"Sequences equals although criteria: {nequals}")
    my_logger.info(md_identical.head())
    
    # Filter fasta
    nrecords = filter_fasta(
        ID_list = ids_to_delete, 
        input_fasta = snakemake.input.fasta, 
        out_fasta = snakemake.output.fasta, grep_invert = True)
    my_logger.info(f"Filtered fasta file: {nrecords} records  - {snakemake.output.fasta} ")

    # Filter metadata
    md = md[~md['NUCCORE_ACC'].isin(ids_to_delete)]
    my_logger.info(f"Filtered metadata file: {len(md.index)} records")
    
    nucc_filt = dfs['nucc'][dfs['nucc'].NUCCORE_ACC.isin(set(md.NUCCORE_ACC))]
    nucc_filt.to_csv(snakemake.output.nucc_filt, index=False)
    my_logger.info(f"# Non-redudant NUCCORE records: {len(nucc_filt.index)}")
    
    # Metadata from identical
    nucc_identical = dfs['nucc'][dfs['nucc'].NUCCORE_ACC.isin(list(md_identical.NUCCORE_ACC))].set_index("NUCCORE_ACC")
    nucc_identical = nucc_identical.reindex(list(md_identical.NUCCORE_ACC))
    nucc_identical['NUCCORE_Identical'] = list(md_identical.NUCCORE_Identical)
    nucc_identical.to_csv(snakemake.output.nucc_identical, index=True)
    my_logger.info(f"# Identical NUCCORE records: {len(nucc_identical.index)}")


