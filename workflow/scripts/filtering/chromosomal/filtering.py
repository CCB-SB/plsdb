
##################################################
# IMPORTS
##################################################
import sys
from os.path  import abspath
sys.path.insert(0, abspath("scripts"))

from utils.utils_my_logger import My_logger
import pandas as pd
import multiprocessing as mp

##################################################
# FUNCTIONS
##################################################

def run_blastn_check(acc, main_fasta, blastn_db, blastn_header, blastn_pident, outfile, my_logger = None):
    """
    Used to run multiple jobs of remote BLASTn searches for chromosome candidates (identified by rMLST analysis).
    The function is called by the pipeline rule "filter3".
    :param acc: Sequence (FASTA) accession
    :param obname: Output (base)name, will be used to create a FASTA and output file for the given accession
    :param main_fasta: FASTA file containing all sequences, will be used to generate a FASTA file containing only the given accession
    :param blastn_cmd: BLASTn CMD (defined in the config file)
    :param blastn_cmd: BLASTn binary (for CMD)
    :param blastn_cmd: BLASTn output header (for CMD, defined in the config file)
    :param blastn_cmd: BLASTn pct. identity (for CMD, defined in the config file)
    """
    import os
    from Bio import SeqIO
    from workflow.scripts.utils.utils import run_cmd

    acc_fasta = f'{outfile}.{acc}.fna'
    acc_ofile = f'{outfile}.{acc}.tsv'

    # Obtain fasta sequence of acc
    with open(main_fasta, 'r') as ifile:
        for record in SeqIO.parse(ifile, 'fasta'):
            if record.id == acc:
                SeqIO.write(record, acc_fasta, 'fasta')
                break
    assert os.path.exists(acc_fasta), f'No file {acc_fasta}'

    cmd = f"""
    blastn -task megablast -query {acc_fasta} -db {blastn_db} \
        -out {acc_ofile} -outfmt '6 {' '.join(blastn_header)}' \
        -perc_identity {blastn_pident} -max_target_seqs 10 -max_hsps 10
    """

    cmd, p_stdout, p_stderr = run_cmd(cmd=cmd, split=True, shell=False)

    return (cmd, p_stdout, p_stderr)
    
def check_blastn_results(acc, acc_fasta, acc_tsv, blastn_header, 
                            blastn_pident, blastn_qcovs, my_logger=None):
    """
    Filter blastn hit results for a given Accession ID (acc)
    
    return: acc if some criteria is met

    """
    import os
    import pandas as pd

    # clean up: rm FASTA
    if os.path.exists(acc_fasta):
        os.remove(acc_fasta)
    # load hits
    try:
        acc_df = pd.read_csv(acc_tsv, sep='\t', header=None, names=blastn_header)
    except pd.errors.EmptyDataError:
        my_logger.info(f'KEEP {acc}: no blastn hits')
        return None
    # filter hits
    acc_df = acc_df.loc[(acc_df['pident'] >= blastn_pident) & (acc_df['qcovs'] >= blastn_qcovs),:]
    # no hits left
    if acc_df.shape[0] == 0:
        my_logger.info(f'KEEP {acc}: no blastn hits left after filtering (pident {blastn_pident}, qcovs {blastn_qcovs})')
        return None
    else:
        top1 = acc_df.iloc[0,:]
        my_logger.info(f"RM {acc}: 1st hit: {top1['sseqid']} [{top1['stitle']}] (evalue {top1['evalue']}, pident {top1['pident']}, qcovs {top1['qcovs']})")
        return acc


##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    from utils.utils import run_cmd, filter_fasta
    
    # Logger
    logger = My_logger(log_filename = snakemake.log[0], logger_name = "rmlst_filtering")
    logger = logger.get_logger()

    # Read plasmid and rmlst dataframes
    pls = pd.read_csv(snakemake.input.pls)
    pls.set_index('NUCCORE_ACC', drop=False, inplace=True, verify_integrity=True)
    logger.info(f'Read in {len(pls.index)} plasmid records\n{pls.head()}')

    rmlst = pd.read_csv(snakemake.input.rmlst, sep='\t', header=0, dtype=str)
    logger.info(f'Read in {len(rmlst.index)} rMLST hits\n{rmlst.head()}')

    #
    ## Count unique loci per query and filter
    #
    pls['rMLST_hits'] = ''
    pls['rMLST_hitscount'] = 0
    for qseqid, df in rmlst[['qseqid', 'slocus']].groupby(['qseqid']):
        hits = sorted(list(set(df['slocus'])))
        pls.loc[qseqid, 'rMLST_hits'] = ';'.join(hits)
        pls.loc[qseqid, 'rMLST_hitscount'] = len(hits)
    pls.to_csv(f'{snakemake.input.pls}.backup', index=False)

    # check those above given cutoff
    check_ids = pls.loc[pls['rMLST_hitscount'] > snakemake.params.cutoff,'NUCCORE_ACC']
    
    #
    ## run BLASTn
    #

    logger.info(f'Running blastn against nuccoredb local using {snakemake.threads} cores')
    args = [(acc, snakemake.input.fasta, snakemake.params.blastndb, 
             snakemake.params.blastn_header, snakemake.params.blastn_pident, 
             snakemake.output.pls, logger
             ) for acc in check_ids]
    
    with mp.Pool(processes = snakemake.threads) as pool:
      try:
          results = pool.starmap(run_blastn_check, args)
      except Exception as e:
          logger.error("Something when wrong with run_blast_check")
          pool.terminate()
          raise e
    logger.info("Done")

    #
    ## check results
    #
    logger.info(f'Checking blastn results using {snakemake.threads} cores')
    args_check = [(acc, f'{snakemake.output.pls}.{acc}.fna' , f'{snakemake.output.pls}.{acc}.tsv',
                  snakemake.params.blastn_header, snakemake.params.blastn_pident, 
                  snakemake.params.blastn_qcovs, logger )for acc in check_ids]
    
    with mp.Pool(processes = snakemake.threads) as pool:
      try:
          drop_indexes = pool.starmap(check_blastn_results, args_check)
      except Exception as e:
          logger.error("Something when wrong with check_blastn_results")
          pool.terminate()
          raise e
    drop_indexes = [acc for acc in drop_indexes if acc]
    pls.drop(index=drop_indexes, inplace=True)
    logger.info(f'Filtered: by number of rMLST and blastn hits: kept {len(pls.index)} records')

    # save metadata
    pls.to_csv(snakemake.output.pls, index=False)

    # Filter fasta
    nrecords = filter_fasta(
        ID_list = list(pls['NUCCORE_ACC']), 
        input_fasta = snakemake.input.fasta,
        out_fasta = snakemake.output.fasta,
        tmp_file = "tmp_ids.txt", grep_invert = False)
    logger.info(f"Filtered fasta file: {int(nrecords)} records  - {snakemake.output.fasta} ")

