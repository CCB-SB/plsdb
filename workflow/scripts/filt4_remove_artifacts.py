import pandas
import argparse
from Bio import SeqIO


##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ifile', '-i', help="Plasmids table", required=True) 
    parser.add_argument('--fasta', '-f', help="Plasmids fasta", required=True)
    parser.add_argument('--ofile', '-o', help="Output file", required=True)
    parser.add_argument('--log', '-l', help="Logger file", required=True)
    return parser


##################################################
# Methods
##################################################

# create dictionary with biosample accession as key and description as additional key
def get_accessions(ifile):
    ptable = pandas.read_csv(ARGS.ifile, sep='\t', header=[0])

    dictionary = dict()
    for index, row in ptable.iterrows():
        isnotna = pandas.notna(row['ACC_BIOSAMPLE'])
        if isnotna:
            key = row['ACC_BIOSAMPLE']
            if key in dictionary.keys():
                dictionary[key].append((row['ACC_NUCCORE'], row['Length_NUCCORE']))
            else:
                dictionary[key] = [(row['ACC_NUCCORE'], row['Length_NUCCORE'])]
        
        key = row['Description_NUCCORE']
        assert pandas.notna(key), "Key for %s is NA" % row['ACC_NUCCORE']
        if key in dictionary.keys():
            dictionary[key].append((row['ACC_NUCCORE'], row['Length_NUCCORE']))
        else:
            dictionary[key] = [(row['ACC_NUCCORE'], row['Length_NUCCORE'])]

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

def get_blastn_results(fpath, output, accs):
    from utils import run_cmd, mkdir
    import os

    cmd_blastn = "blastn -query {query} -subject {subject} -outfmt '6 qseqid sseqid qlen slen length pident qstart qend evalue' -max_hsps 5 >> {out}"
    idpath = '../results/data/tmp/' 
    mkdir(idpath)
    # split all-fastas file into seperate files
    for record in SeqIO.parse(fpath, 'fasta'):
        #idpath = os.path.abspath('/tmp/' + record.id + '.fasta')
        f = record.id + '.fasta'
        
        #print(idpath)
        idfile = open(idpath + f, 'w')
        idfile.write(">"+str(record.id)+'\n'+str(record.seq))
        idfile.close()

    # step through all accessions grouped together in previous step
    for samples in accs:
        for i in range(0, len(samples)):
            for j in range(i+1, len(samples)):
                query = '../results/data/tmp/' + str(samples[j]) + '.fasta'
                subject = '../results/data/tmp/' + str(samples[i]) + '.fasta'
                # run blastn between query and subject and append output file
                cmd = cmd_blastn.format(subject=subject, query=query, out=output)
                cmd, cmd_s, cmd_o = run_cmd(cmd)
                assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd, cmd_s, cmd_o)


def find_artifacts(blastn_results, ptable, output):
    # open blastn file created in previous step
    blastn = pandas.read_csv(blastn_results, sep='\t', names=['qseqid', 'sseqid', 'qlen', 'slen', 'length', 'pident', 'qstart', 'qend', 'evalue'])

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

        value = length/qlen * pident

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

    plasmids = pandas.read_csv(ptable, sep='\t', header=[0])

    # keep non-artifacts in table
    pls = plasmids.loc[~plasmids['ACC_NUCCORE'].isin(artifacts)]
    pls['inclusions'] = ''
    pls.set_index('ACC_NUCCORE', inplace=True)
    for i, inc in enumerate(incl):
        pls.loc[inc,'inclusions'] += artifacts[i] + ';'
        
    pls.reset_index(inplace=True)
    pls.set_index('UID_NUCCORE', inplace=True)
    pls.to_csv(output, sep='\t', index=True)
    return(artifacts)


##################################################
# MAIN
##################################################
if __name__ == "__main__":
    from utils import run_cmd, setup_logger
    import logging
    # Args
    ARGS = get_arg_parser().parse_args()
    
    logger = setup_logger(logging.INFO, ARGS.log)

    accessions = get_accessions(ARGS.ifile)
    logger.info("Found {} groups with same NUCCORE description or BIOSAMPLE ID".format(len(accessions)))

    logger.info("Start blastn to compare plasmids in same group")
    blastn_results = 'tmp.blastn.tsv'
    get_blastn_results(ARGS.fasta, blastn_results, accessions)
    
    logger.info("Process blastn results to identify artifacts")
    artifacts = find_artifacts(blastn_results, ARGS.ifile, ARGS.ofile)
    logger.info("Found {} artifacts: {}".format(len(artifacts), artifacts))
    # remove dowloaded fastas and tmp.blastn.tsv
    cmd = "rm -r ../results/data/tmp & rm tmp.blastn.tsv"
    cmd, cmd_s, cmd_o = run_cmd(cmd)
    assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd, cmd_s, cmd_o)
    logger.info("Temporary files were removed")
    




    
