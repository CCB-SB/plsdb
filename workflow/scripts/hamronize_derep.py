import pandas as pd
import numpy as np
from utilsmeta.my_logger import My_logger

# Logging
logger = My_logger(log_filename = snakemake.log[0], logger_name = "process_calculate_gc")
my_logger = logger.get_logger()

df = pd.read_table(snakemake.input[0], low_memory=False)
my_logger.info(f"# Initial rows: {len(df.index)}")

my_logger.info(f"Priority: AMRFinderPlus, RGI")

# Detect equals gene symbols. Example: OXA-48 (rgi) and blaOXA-48 (ncbi)
ncbi_df = df[df['analysis_software_name'] == 'amrfinderplus']
bla_ncbi_genes = set(i.replace('bla', '') for i in set(ncbi_df['gene_symbol']) if i.startswith('bla'))
df['gene_symbol2'] = np.where((df['analysis_software_name'] == 'rgi') & (df['gene_symbol'].isin(bla_ncbi_genes)),
         'bla'+df['gene_symbol'].str.lower(), df['gene_symbol'].str.lower())

# First filtering
cols_dup = ['input_sequence_id', 'gene_symbol2', 'strand_orientation',
            'input_gene_start', 'input_gene_stop']
df = df.sort_values(by='analysis_software_name', ascending=True) # amrfinder, rgi priority list
df.drop_duplicates(subset=cols_dup, keep='first', inplace=True)
df.reset_index(drop=True, inplace=True)
my_logger.info(f"# 1st Deduplicated rows: {len(df.index)}")

# Aplying gene start and stop to range +-10
range_nt = snakemake.params.range_nt
my_logger.info(f"Applying range in input_gene_start and input_gene_stop: +-{range_nt}")
dup = df[df.duplicated(subset=cols_dup[:-2], keep=False)]

## TOOO SLOWWWWW
my_logger.info(f"# Detected possible duplicates: {len(dup.index)}")

# Compute distances
index_to_del = []
for seq in set(dup['input_sequence_id']):
    m = {}
    for i in ['start', 'stop']:
        # Compute pairwise distances
        dup_sub = dup[dup['input_sequence_id']==seq]
        # my_logger.debug(dup_sub[['input_sequence_id', 'gene_symbol2', 'strand_orientation',
            # 'input_gene_start', 'input_gene_stop']])
        matrix = dup_sub[f'input_gene_{i}'].values.reshape(-1, 1) - dup_sub[f'input_gene_{i}'].values
        # my_logger.debug(matrix)
        
        # Set nan upper triangle and diagonal to remove redundancy
        triu = np.tri(matrix.shape[0], matrix.shape[1], k=-1) # fill with 0
        triu[triu==0] = np.nan
        matrix = matrix * triu
        # my_logger.debug(matrix)
        
        # Set nan distances above threshold
        matrix = np.where(abs(matrix) > range_nt, np.nan, matrix)
        m[i] = matrix
        # my_logger.debug(m[i])

    # Start and Stop distance must not exceed threshold
    sum_matrix = abs(m['start']) + abs(m['stop'])
    # my_logger.debug(sum_matrix)
    sum_matrix = np.where(sum_matrix > range_nt*2, np.nan, sum_matrix)
    # my_logger.debug(sum_matrix)

    # Get dups
    dups = np.transpose(np.where(~np.isnan(sum_matrix)))
    # my_logger.info(f"# Duplicated pairs: {len(dups)}")
    # my_logger.info(dups)

    # Get the highest index because corresponds to rgi
    seq_position_to_del  = [max(pair) for pair in dups]
    index_to_del_sub = list(dup_sub.iloc[seq_position_to_del].index)
    # my_logger.debug(index_to_del_sub)
    index_to_del.extend(index_to_del_sub)

df = df.loc[df.index.delete(index_to_del), :]
df.to_csv(snakemake.output[0], index=False, sep='\t')
my_logger.info(f"# 2nd Deduplicated rows: {len(df.index)}")


