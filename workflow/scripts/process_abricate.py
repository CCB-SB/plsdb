# NOTE: Using a procedure analoguous to that of PlasmidFinder (https://bitbucket.org/genomicepidemiology/plasmidfinder)
#   1) Call Blaster (from CGE core module): call BLAST and pre-process hits (best hit per subject)
#   2) Filter results by identity and coverage; from overlapping hits keep only the best hit
#   Differences/modifications: see tag "CHANGED"


import pandas as pd
import os
from Bio import SeqIO
from utils_my_logger import My_logger
import multiprocessing as mp



##################################################
# ARGS
##################################################
def get_arg_parser():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--fna', '-f', help="Fna file.", required=True)
    parser.add_argument('--abr', help="File for updating abricate.")
    parser.add_argument('--ofile', '-o', help="Output file.", required=True)
    parser.add_argument('--cov', help="Abricate cov parameter.", type=int ,required=True)
    parser.add_argument('--ident', help="Abricate ident parameter.", type=int, required=True)
    parser.add_argument('--env_path', help="Path to conda environment.", required=True)
    parser.add_argument('--wildcard', help="Wildcard to be used", required=True)
    parser.add_argument('--log', help="Log File", required=False)
    parser.add_argument("--cores", type = int, default = 1, required = False,
                        help="Cores to use" )
    return parser

##################################################
# FUNCTIONS
##################################################

def call_blaster(record, db, db_path,  cov, ident, my_logger = None):
    import os
    from Bio import SeqIO
    from blaster import Blaster

    tmp_fasta   = f'{record.id}.{db}.fasta.tmp'
    tmp_blaster = f'{record.id}.{db}.blaster.tmp'

    # Create FASTA
    SeqIO.write(record, tmp_fasta, 'fasta')

    # Call Blaster
    if os.path.exists(tmp_blaster):
        os.remove(tmp_blaster)
    try:
        results = Blaster(
            inputfile=tmp_fasta,
            databases=['sequences'], # created by Abricate, *.fsa link
            db_path=db_path,
            out_path=tmp_blaster,
            min_cov=cov/100.0, # expects 0.0 .. 1.0
            threshold=ident/100.0, # expects 0.0 .. 1.0
            cut_off=False
        )
    except Exception as e:
         my_logger.error(f"Blaster cmd: inputfile={tmp_fasta}")
         raise e

    # Collect results
    if results.results['sequences'] == 'No hit found':
        os.remove(tmp_fasta)
        os.remove(tmp_blaster)
        return

    df = []
    for result in results.results['sequences'].values():
        df.append({
            'qseqid': result['contig_name'].split(' ')[0],
            'qstart': result['query_start'],
            'qend': result['query_end'],
            # 'qlen':,
            'sseqid': result['sbjct_header'].split(' ')[0],
            'stitle': result['sbjct_header'],
            'sstart': result['sbjct_start'],
            'send':result['sbjct_end'],
            'slen': result['sbjct_length'],
            'evalue': result['evalue'],
            'pident': result['perc_ident'],
            'length': result['HSP_length'],
            'gaps': result['gaps'],
            'cov': result['perc_coverage']
        })

    # rm tmp files
    os.remove(tmp_fasta)
    os.remove(tmp_blaster)
    
    return df

def process_blaster(result, cov, pident):
    import pandas as pd
    import copy

    # no results
    if not result:
        return None, None

    # convert to table
    result = pd.DataFrame(result)

    # filter by coverage and identity
    result = result.loc[(result['cov'] >= cov) & (result['pident'] >= pident),:]

    # no results left
    if len(result.index) == 0:
        return None, None

    # compute score from identity and coverage
    result = result.assign(ic_score = result['pident'] * result['cov'])

    # sort by query start: small -> large (ascending)
    result.sort_values(by=['qstart'], ascending=[True], inplace=True)
    result_all = copy.copy(result)

    # overlapping hits
    keep_hit = [True] * result.shape[0] # True = keep, False = remove
    current_i   = 0
    current_end = result['qend'].values[0]
    current_ics = result['ic_score'].values[0]
    for i in range(0, result.shape[0]-1):
        next_start = result['qstart'].values[i+1]
        next_end   = result['qend'].values[i+1]
        next_ics   = result['ic_score'].values[i+1]

        if next_start <= current_end: # overlapping hits; CHANGED: <=, NOT <
            # remove hit with lower identity-coverage score
            if current_ics < next_ics:
                keep_hit[current_i] = False
                # update current
                current_i   = i+1
                current_end = next_end
                current_ics = next_ics
            else:
                    keep_hit[i+1] = False
        else: # no overlap, update current; CHANGED: otherwise if no overlap always comparing to the latest "current"
                current_i   = i+1
                current_end = next_end
                current_ics = next_ics
    result = result.loc[keep_hit,:]

    # rm coverage-identity score
    result = result.drop(['ic_score'], axis=1, inplace=False) # inplace=True produces SettingWithCopyWarning

    return result, result_all

def blaster_abricate(record, db, db_path,  cov, ident, my_logger = None):
    df = call_blaster(record = record, db = db, db_path = db_path,
                           cov = cov, ident=ident, my_logger=my_logger)
    df, df_all = process_blaster(result = df, cov = cov, pident=ident)

    return df, df_all
##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    ARGS = get_arg_parser().parse_args()
    
    # Logger
    logger = My_logger(log_filename = ARGS.log, logger_name = "process_abricate")
    logger = logger.get_logger()

    db_path = os.path.join(ARGS.env_path, 'db', ARGS.wildcard)
    logger.info(f'Blaster database in {db_path}')

    # call Blaster to get and pre-process BLAST hits
    logger.info(f"Call Blaster for each record in {ARGS.fna}")
    args = [(record, ARGS.wildcard, db_path,
             ARGS.cov, ARGS.ident, logger) for record in SeqIO.parse(ARGS.fna, 'fasta')]
    
    with mp.Pool(processes = ARGS.cores) as pool:
        try:
            results = pool.starmap(blaster_abricate, args)
        except Exception as e:
            logger.error("Something when wrong with blaster_abricate")
            pool.terminate()
            raise e
    
    df_all = [df_all_sub for df_sub, df_all_sub in results if df_all_sub is not None]
    df = [df_sub for df_sub, df_all_sub in results if df_sub is not None]
    
    if not df_all and not df:
        logger.info("No results, creating empty dataframes...")
        df = pd.DataFrame(columns=["qseqid","qstart","qend",
                                   "sseqid","stitle","sstart",
                                   "send","slen","evalue","pident",
                                   "length","gaps","cov","sseqdb"])
        df_all = pd.DataFrame(columns=["qseqid","qstart","qend",
                                   "sseqid","stitle","sstart",
                                   "send","slen","evalue","pident",
                                   "length","gaps","cov", "ic_core",
                                   "sseqdb"])
        df_all.to_csv(f"{ARGS.ofile}.all", index=False)
        df.to_csv(ARGS.ofile, index=False)
        exit()
    
    df_all = pd.concat(df_all, axis=0, sort=False)
    df = pd.concat(df, axis=0, sort=False)


    # stitle -> DB + ID
    # title format: <BLAST DB ID> <DB>~~~<name>~~~<ID> <name>
    # example: gnl|BL_ORD_ID|578 argannot~~~(Bla)blaFRI-3~~~KY524440:1-885 (Bla)blaFRI-3
    df_all['sseqdb'] = df_all['stitle'].map(lambda x: x.split(' ')[1].split('~~~')[0])
    df_all['sseqid'] = df_all['stitle'].map(lambda x: ', '.join(x.split(' ')[1].split('~~~')[1:]))
    df['sseqdb']     = df['stitle'].map(lambda x: x.split(' ')[1].split('~~~')[0])
    df['sseqid']     = df['stitle'].map(lambda x: ', '.join(x.split(' ')[1].split('~~~')[1:]))
    
    # save
    df_all.to_csv(f"{ARGS.ofile}.all", index=False)
    df.to_csv(ARGS.ofile, index=False)
