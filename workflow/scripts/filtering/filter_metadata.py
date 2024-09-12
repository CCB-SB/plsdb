#!/usr/bin/env python
# coding: utf-8


##################################################
# FUNCTIONS
##################################################


def filter_by_metadata(pls, filter_regex, my_logger=None):
    """
        Filter plasmid based on metadata attributes:
            - Duplicated entry (NUCCORE_DuplicatedEntry)
            - Description (regex)
            - Assembly Completeness (Complete)
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

    return pls


##################################################
# MAIN
##################################################

if __name__ == "__main__":
    from utilsmeta.my_logger import My_logger
    from utilsmeta.utils import merge_NABT, read_NABT
    import pandas as pd

    # Logging
    logger = My_logger(log_filename = snakemake.log[0], logger_name = "plasmid_metadata_filtering")
    my_logger = logger.get_logger()

    dfs = read_NABT(nucc_path=snakemake.input.nucc, 
              bio_path=snakemake.input.bio, 
              ass_path=snakemake.input.ass,
              tax_path=snakemake.input.tax)


    # Merge Nuccore + Assembly
    nucc_ass_df = merge_NABT(nucc=dfs['nucc'], 
                             ass=dfs['ass'])
    print(nucc_ass_df.head())
    
    # Filter by metadata
    pls = filter_by_metadata(pls = nucc_ass_df, 
                             filter_regex = snakemake.params.filter_regex,
                    my_logger = my_logger)
    
    # Filter
    nucc_df = dfs['nucc'][dfs['nucc'].loc[:, 'NUCCORE_UID'].isin(list(pls['NUCCORE_UID']))]
    ass_df = dfs['ass'][dfs['ass'].loc[:, 'ASSEMBLY_UID'].isin(set(pls['ASSEMBLY_UID']))]
    bio_df = dfs['bio'][dfs['bio'].loc[:, 'BIOSAMPLE_UID'].isin(set(pls['BIOSAMPLE_UID']))]
    tax_df = dfs['tax'][dfs['tax'].loc[:, 'TAXONOMY_UID'].isin(set(pls['TAXONOMY_UID']))]

    nucc_df.drop(columns='Unnamed: 0').to_csv(snakemake.output.nucc, index=False)
    ass_df.drop(columns='Unnamed: 0').to_csv(snakemake.output.ass, index=False)
    bio_df.drop(columns='Unnamed: 0').to_csv(snakemake.output.bio, index=False)
    tax_df.drop(columns='Unnamed: 0').to_csv(snakemake.output.tax, index=False)

