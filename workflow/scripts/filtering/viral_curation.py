##################################################
# IMPORTS
##################################################


from utilsmeta.my_logger import My_logger

suspicious_genes = [
    "carbohydrate kinase", "carbohydrate-kinase", "glycosyltransferase",
    "glycosyl transferase", "glycosyl transferaseendonuclease",
    "nucleotide sugar epimerase", "nucleotide sugar-epimerase",
    "nucleotide-sugar epimerase", "nucleotide-sugar-epimerase",
    "nucleotidyltransferase", "nucleotidyl transferase",
    "nucleotidyl-transferase","plasmid stability",
    "endonuclease"
]

##################################################
# MAIN
##################################################
if __name__ == "__main__":
    
    import numpy as np
    import pandas as pd
    
    # Logging
    logger = My_logger(log_filename = snakemake.log[0], 
                       logger_name = "viral_filtering")
    my_logger = logger.get_logger()

    #
    ## Load df
    #

    # Extract sequence identifier
    vir = pd.read_csv(snakemake.input.virsorter2, sep='\t', header=0)
    my_logger.info(f"# Virsorter2 results = {len(vir.index)}")

    checkv = pd.read_csv(snakemake.input.checkv, sep="\t", header = 0)
    checkv = checkv.rename(columns={'contig_id': 'seqname'})
    my_logger.info(f"# Checkv results = {len(checkv.index)}")

    dramv = pd.read_csv(snakemake.input.dramv, sep='\t', header=0)
    seqids = [i[0] for i in dramv['gene'].str.split('__')]
    dramv['seq_id'] = seqids
    my_logger.info(f"# DRAMv results = {len(dramv.index)}")

    # Merge 

    df = vir.merge(checkv, how='outer', on='seqname',
                   validate="1:1").drop_duplicates().reset_index(drop=True)
    seqids = [i[0] for i in df['seqname'].str.split('|')]
    df['seq_id'] = seqids

    #
    ## Filter
    #

    my_logger.info(f"""Now we will evaluate the results following the sequential criteria:
                    - likely_viral_1 : >0 viral_genes by Checkv
                    - likely_viral_2: 0 viral_genes & (0 host_genes by Checkv or >2 hallmarks by Virsorter2 or >0.95 max_score by Virsorter)
                    - mostly_non-viral: 0 viral_genes & >=2 host_genes
                    - mostly_mge: 0 viral_genes & 1 host_genes
                    - undefined: other cases
                   
                   More info: https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3
                   """)

    # >0 viral_gene = VIRAL
    df['eval_1'] = np.where(df['viral_genes'] > 0, "likely_viral_1",
                np.where((df['viral_genes'] == 0) & ((df['host_genes'] == 0) | (df['hallmark'] > 2) | (df['max_score'] > 0.95) ), "likely_viral_2",
                         np.where((df['viral_genes'] == 0) & (df['host_genes'] >= 2), "mostly_non-viral",
                                  np.where((df['viral_genes']==0) & (df['host_genes'] == 1), "mostly_mge",
                                    "undefined")
                                  )
                )
    )
    my_logger.info(f"Viral evaluation: {df.value_counts(subset='eval_1')}")

    # Scan putative viral in DRAMv results for false positives
    putative_viral = set(df[df['eval_1'] == 'likely_viral_2']['seq_id'])
    dramv_sub = dramv[dramv['seq_id'].isin(putative_viral)]

    false_positives = list(dramv_sub[dramv_sub['gene_description'].str.lower(
        ).isin(suspicious_genes) | dramv_sub['module'].str.lower(
          ).isin(suspicious_genes)]['seq_id'])
    
    df = df[~df['seq_id'].isin(false_positives)]
    my_logger.info("Scanning for possible false_positives in the likely_viral_2 category using DRAMv results")
    my_logger.info(f"# False positives: {len(false_positives)}")

    # Relabel
    df['putative_viral'] = df['eval_1']
    df['putative_viral'] = np.where(df['putative_viral'].str.contains('likely_viral'), True, False)
    df = df.rename(columns={'eval_1': 'VIRAL_category',
               'putative_viral': 'VIRAL_putative',
               'seqname': "VIRAL_contig",
               'seq_id': 'NUCCORE_ACC'})
    my_logger.info(f"Final putative_viral annotation, merging categories: {df.value_counts(subset='VIRAL_putative')}")
    my_logger.info(df.head())
    print(df.columns)
    
    # Locate which are different
    d = {}
    for acc, cat in zip(df['NUCCORE_ACC'], df['VIRAL_category']):
      d.setdefault(acc, {}).setdefault('VIRAL_category', set()).add(cat)
    
    for key, value in d.items():
      if len(value['VIRAL_category']) > 1:
        if 'likely_viral' in list(value['VIRAL_category'])[0] and list(value['VIRAL_category'])[1]:
           d[key]['VIRAL_category'] = 'likely_viral_1-2'
        else:
           d[key]['VIRAL_category'] = 'undefined'
      else:
        d[key]['VIRAL_category'] = list(value['VIRAL_category'])[0]
    
    df_nucc = pd.DataFrame.from_dict(d, orient='index').reset_index(names="NUCCORE_ACC")
    df_nucc['VIRAL_putative'] = np.where(df_nucc['VIRAL_category'].str.contains('likely_viral'), True, False)
    print(df_nucc)
  
    df_nucc.to_csv(snakemake.output.pls, index=False)
    df.to_csv(snakemake.output.viral_tab, index=False)