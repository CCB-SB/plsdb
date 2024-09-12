from utilsmeta.my_logger import My_logger
from Bio import SeqIO



##################################################
# Methods
##################################################

# create dictionary with biosample accession as key and description as additional key
def get_accessions(ifile):
    import pandas as pd

    ptable = pd.read_csv(ifile)

    dictionary = dict()
    for index, row in ptable.iterrows():
        isnotna = pd.notna(row['BIOSAMPLE_ACC'])
        if isnotna:
            key = row['BIOSAMPLE_ACC']
            if key in dictionary.keys():
                dictionary[key].append((row['NUCCORE_ACC'], row['NUCCORE_Length']))
            else:
                dictionary[key] = [(row['NUCCORE_ACC'], row['NUCCORE_Length'])]
        
        key = row['NUCCORE_Description']
        assert pd.notna(key), f"Key for {row['NUCCORE_ACC']} is NA"
        if key in dictionary.keys():
            dictionary[key].append((row['NUCCORE_ACC'], row['NUCCORE_Length']))
        else:
            dictionary[key] = [(row['NUCCORE_ACC'], row['NUCCORE_Length'])]

    dictionary = {key: item for key, item in dictionary.items() if len(item) > 1}

    accessions = []
    for key, item in dictionary.items():
        # sort each array in dict according to plasmid length; largest first
        item.sort(key=lambda x:-x[1])
        # store ids of samples with same key (biosample accession/description) together
        sample = []
        for i in item:
            sample.append(i[0])
        accessions.append(sample)

    return accessions

def get_blastn_results(fpath, output, accs, tmp_dir, cores = 10):
    from workflow.scripts.utils.utils import run_cmd, mkdir
    from os.path import join, exists
    from os import remove
    import multiprocessing as mp

    # split all-fastas file into seperate files
    mkdir(tmp_dir)
    for record in SeqIO.parse(fpath, 'fasta'):
        f = f'{record.id}.fasta'
        
        with open(join(tmp_dir, f), 'w') as idfile:
            idfile.write(">"+str(record.id)+'\n'+str(record.seq))

    # step through all accessions grouped together in previous step
    args = []
    tmp_prefix = "tmp_blastn0"
    for samples in accs:
        for i in range(0, len(samples)):
            for j in range(i+1, len(samples)):
                query = join(tmp_dir, f"{str(samples[j])}.fasta")
                subject = join(tmp_dir, f"{str(samples[i])}.fasta")
                output_tmp = f"{tmp_prefix}{i}"
                args.append((query, subject, output_tmp))
  
    # run blastn between query and subject and append output file
    cmd_blastn = [f"blastn -query {query} -subject {subject} -outfmt '6 qseqid sseqid qlen slen length pident qstart qend evalue' -max_hsps 5 > {output_tmp}" for query, subject, output_tmp in args]
    cmd_args = [(cmd, False, True) for cmd in cmd_blastn]
    
    with mp.Pool(processes = cores) as pool:
        try:
            results = pool.starmap(run_cmd, cmd_args)
        except Exception as e:
            logger.error("Something when wrong with get_blastn_results")
            pool.terminate()
            raise e
    
    run_cmd(f"cat {tmp_prefix}* >> {output}", split=False, shell=True)
    [remove(out_tmp) for query, subject, out_tmp in args if exists(out_tmp)]
 


def find_artifacts(blastn_results, ptable, output, my_logger = None):
    import pandas as pd

    # open blastn file created in previous step
    blastn = pd.read_csv(blastn_results, 
                         sep='\t', 
                         names=['qseqid', 'sseqid', 'qlen', 
                                'slen', 'length', 'pident', 
                                'qstart', 'qend', 'evalue'])

    artifacts = [] 
    d = dict()
    # fill dictionary with query-subject keys
    for index,row in blastn.iterrows():
        key = (row['qseqid'],row['sseqid'])
        item = (row['qstart'], row['qend'], row['qlen'], row['length'], row['pident'])
        if key in d:
            d[key].append(item)
        else:
            d[key] = [item]

    tophits = dict()
    for key, items in d.items():
        # consider that plasmids are circular and blastn handles only linear alignments, 2 HSP can form top alignment
        items = sorted(items)
        start = 0
        end = 0
        length = 0
        pident = 0
        qlen = 0
        for item in items:
            if item[0] == end + 1:
                end = item[1]
                qlen = item[2]
                length = end
                if pident == 0:
                    pident = item[4]
                else:
                    pident = (pident + item[4])/2
            elif length < item[3]:
                start = item[0]
                end = item[1]
                qlen = item[2]
                length = item[3]
                pident = item[4]
        try:
            value = length/qlen * pident
        except Exception as e:
            my_logger.error(f"blastn df:\n{blastn.head()}")
            my_logger.error(f"qlen = {qlen}")
            raise e

        # search for each plasmid the highest coverage value
        if key[0] in tophits:
            if value > tophits[key[0]][1]:
                tophits[key[0]] = (key[1], value, qlen)
        else:
            tophits[key[0]] = (key[1], value, qlen)

    # store artifacts
    incl = []
    for key, item in tophits.items():
        if item[1] == 100:
            artifacts.append(key)
            incl.append(item[0])

    plasmids = pd.read_csv(ptable)
    my_logger.info(f"Plasmid records before artifact filtering : {len(plasmids.index)}")

    # keep non-artifacts in table
    pls = plasmids.loc[~plasmids['NUCCORE_ACC'].isin(artifacts)]
    pls['inclusions'] = ''
    pls.set_index('NUCCORE_ACC', inplace=True)
    for i, inc in enumerate(incl):
        pls.loc[inc,'inclusions'] += artifacts[i] + ';'
        
    pls.reset_index(inplace=True)
    pls.set_index('NUCCORE_UID', inplace=True)
    pls.to_csv(output, index=True)

    my_logger.info(f"Plasmid records after artifact filtering : {len(pls.index)}")

    return(artifacts)


##################################################
# MAIN
##################################################
if __name__ == "__main__":
    from workflow.scripts.utils.utils import run_cmd, filter_fasta
    import shutil
    import sys

    # Args
    ARGS = get_arg_parser().parse_args()
    
    # Logger
    logger = My_logger(log_filename = ARGS.log, logger_name = "remove_artifacts")
    logger = logger.get_logger()

    accessions = get_accessions(ARGS.ifile)
    logger.info(f"Found {len(accessions)} groups with same NUCCORE description or BIOSAMPLE ID")
    
    if len(accessions)==0:
        shutil.copy2(ARGS.ifile, ARGS.out_md)
        sys.exit(0)

    logger.info(f"Start blastn to compare plasmids in same group using {ARGS.cores} cores")
    tmp_dir = "../results/tmp"
    blastn_results = 'tmp.blastn.tsv'
    get_blastn_results(ARGS.fasta, blastn_results, accessions, tmp_dir=tmp_dir,
                       cores = ARGS.cores)
    
    # Find artifacts and filter plasmid metadata
    logger.info("Process blastn results to identify artifacts")
    artifacts = find_artifacts(blastn_results, ARGS.ifile, ARGS.out_md,
                               my_logger=logger)
    logger.info(f"Found {len(artifacts)} artifacts: {artifacts}")

    # Filter fasta
    if artifacts:
        nrecords = filter_fasta(
            ID_list = artifacts, 
            input_fasta = ARGS.fasta, out_fasta = ARGS.out_fasta,
            tmp_file = "tmp_ids.txt", grep_invert = True)
        logger.info(f"Filtered fasta file: {int(nrecords)} records  - {ARGS.out_fasta} ")
    else:
        shutil.copy2(ARGS.fasta, ARGS.out_fasta)
        
    # remove downloaded fastas and tmp.blastn.tsv
    cmd = f"rm -r {tmp_dir} & rm {blastn_results}"
    run_cmd(cmd, split=False, shell=True)
    logger.info("Temporary files were removed")
