import pandas as pd
from utils_my_logger import My_logger

##################################################
# ARGS
##################################################
def get_arg_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pls', help="Plasmids DIS ONT", required=True)
    parser.add_argument('--emb', help="Master MASH UMAP.", required=True)
    parser.add_argument('--abr', help="Master Abricate.", required=True)
    parser.add_argument('--pmlst', help="Master pMLST.", required=True)
    parser.add_argument('--ofile', '-o', help="Output file.", required=True)
    parser.add_argument('--log', help="Log File", required=False)
    return parser


##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    ARGS = get_arg_parser().parse_args()
    # Logger
    logger = My_logger(log_filename = ARGS.log, logger_name = "process_infotable")
    logger = logger.get_logger()

    # plasmids
    pls = pd.read_csv(ARGS.pls, dtype=str)
    pls.set_index(keys='NUCCORE_ACC', drop=False, inplace=True, verify_integrity=True)
    logger.info(f'Read in {len(pls.index)} plasmid records\n{pls.head()}')

    # embedding: NOTE
    emb = pd.read_csv(ARGS.emb, sep='\t', header=0, index_col='ID', dtype=str)
    logger.info(f'Read in embedding for {len(emb.index)} records\n{emb.head()}')

    # PlasmidFinder
    abr = pd.read_csv(ARGS.abr)
    abr = abr.loc[abr['sseqdb'] == 'plasmidfinder',['qseqid', 'sseqid']]
    logger.info(f'Read in {len(abr.index)} PlasmidFinder hits\n{abr.head()}')
    pf = []; pf_index = []
    for aggr_id, aggr_df in abr.groupby(['qseqid']):
        assert len(set(aggr_df['sseqid'])) == len(aggr_df['sseqid'])
        pf_index.append(list(aggr_df['qseqid'])[0])
        pf.append({'plasmidfinder': '|'.join(aggr_df['sseqid'])})
    pf = pd.DataFrame(pf, index=pf_index)
    logger.info(f'Aggr PlasmidFinder hits\n{pf.head()}')

    # pMLST
    pmlst_df = pd.read_csv(ARGS.pmlst, header=0)
    logger.info(f'Read in {len(pmlst_df.index)} pMLST records\n{pmlst_df.head()}')
    pmlst = []; pmlst_index = []
    for aggr_id, aggr_df in pmlst_df.groupby(['ID']):
        assert len(set(aggr_df['pmlst'])) == len(aggr_df['pmlst'])
        pmlst_index.append(list(aggr_df['ID'])[0])
        pmlst.append({'pmlst': '|'.join(aggr_df['pmlst'])})
    pmlst = pd.DataFrame(pmlst, index=pmlst_index)
    logger.info(f'Aggr pMLST hits\n{pmlst.head()}')

    # process source strings
    source_map = {'insd': 'INSDC', 'refseq': 'RefSeq'}
    pls['NUCCORE_Source'] = pls['NUCCORE_Source'].map(lambda x: source_map[x])

    # Merge data
    pls = pd.merge(left=pls, right=emb, how='left',
        left_index=True, right_index=True, sort=False)
    logger.info('Added embedding')

    pls = pd.merge(left=pls, right=pf, how='left',
        left_index=True, right_index=True, sort=False)
    logger.info('Added PlasmidFinder hits')

    pls = pd.merge(left=pls, right=pmlst, how='left',
        left_index=True, right_index=True, sort=False)
    logger.info('Added pMLST hits')

    # save
    logger.info(f'Merging completed. Total records: {len(pls.index)}')
    pls.to_csv(ARGS.ofile, sep='\t', header=True,
        index=False, index_label=False)