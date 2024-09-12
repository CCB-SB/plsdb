# IMPORTS
##################################################
import sys
from os.path  import abspath
sys.path.insert(0, abspath("scripts"))

##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    import pandas as pd
    from utils.utils_my_logger import My_logger
    from utils.utils import run_cmd

    # Logger
    logger = My_logger(log_filename = snakemake.log[0], logger_name = "rmlst_blastn")
    logger = logger.get_logger()

    # blastn search
    logger.info(f"Searching plasmid sequences against rmlst database using Blastn...")
    cmd = f"""blastn -query {snakemake.input.fasta} -db {snakemake.input.blastdb} \
        -task blastn -perc_identity {snakemake.params.ident}\
        -out {snakemake.output.tsv} \
        -outfmt '6 {' '.join(snakemake.params.header)}' \
        -num_threads {snakemake.threads}\
        -evalue 1E-20 -culling_limit 1"""
    logger.info(cmd)
    run_cmd(cmd, split=True, shell=False)
    logger.info("Blastn completed")

    # read in results
    df = pd.read_csv(snakemake.output.tsv, sep='\t', header=None, 
                     names=snakemake.params.header)
    logger.info(f"# Blast hits = {len(df.index)}")

    # compute coverage
    df['cov'] =  100 * (df['length'] - df['gaps']) / df['slen']
    
    # filter by coverage
    df = df.loc[df['cov'] >= snakemake.params.cov,:]
    logger.info(f"# Blast hits with coverage > {snakemake.params.cov} = {len(df.index)}")


    # rMLST locus
    df['slocus'] = df['sseqid'].map(lambda x: x.split('_')[0])

    # save
    df.to_csv(
        path_or_buf=snakemake.output.tsv,
        sep='\t',
        header=True,
        index=False,
        index_label=False
    )

    logger.info(f"File generated: {snakemake.output.tsv}")
