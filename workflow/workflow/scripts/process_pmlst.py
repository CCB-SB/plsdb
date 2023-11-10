
import pandas
import os
import re
from Bio import SeqIO
from utils import run_cmd
from utils_my_logger import My_logger

##################################################
# ARGS
##################################################
def get_arg_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', help="pLMST database", required=True)
    parser.add_argument('--fasta', help="Fasta file.", required=True)
    parser.add_argument('--pf', help="Plasmid finder.", required=True)
    parser.add_argument('--ofile', '-o', help="Output file.", required=True)
    parser.add_argument('--map', help="pLMST map from config file.", nargs='+', required=True)
    parser.add_argument('--cmd', help="pMLST command from config file.", required=True)
    parser.add_argument('--env_path', help="Path to conda environment.", required=True)
    parser.add_argument('--log', help="Log File", required=False)
    return parser

##################################################
# FUNCTIONS
##################################################

def reproc_mlst_scheme_name(scheme_name):
    """
    Reverts the changes made by proc_mlst_scheme_name()
    """
    import re
    return re.sub('_', '/', re.sub('__', ' ', scheme_name))

def process_pmlst_hits(f, pmlst_db_path, logger):
    """
    Parse hits from mlst tool (https://github.com/tseemann/mlst)
    :param f: results file
    :param pmlst_db_path: path to the dir containing the pMLST files
    :param logger: logger to be used
    :return list of hits, hit = "<scheme name>(ST): <allele hits>"
    """
    import os
    hits = []
    profiles = {}
    with open(f, 'r') as ifile:
        for line in ifile:
            line = line.rstrip('\n').split('\t')
            # no hit
            if line[1] == "-":
                continue
            # process hit
            if line[2] != "-":
                pfile = os.path.join(pmlst_db_path, line[1], line[1] + '.txt')
                if line[1] not in profiles:
                    if os.path.exists(pfile + '.dummy'):
                        profiles[line[1]] = 'dummy'
                    elif os.path.exists(pfile + '.old'):
                        profiles[line[1]] = pandas.read_csv(pfile + '.old', sep='\t', header=0)
                    else:
                        profiles[line[1]] = None
                if profiles[line[1]] is None:
                    pass
                elif isinstance(profiles[line[1]], str) and profiles[line[1]] == 'dummy': # dummy STs
                    logger.info('Dummy ST: {}'.format(line[1]))
                    line[2] = '-'
                elif 'oldST' in list(profiles[line[1]]): # ST mapping
                    logger.info('Mapped ST: {}'.format(line[1]))
                    line[2] = (profiles[line[1]].loc[profiles[line[1]]['ST'] == int(line[2]),'oldST'].tolist())[0]
            line[1] = reproc_mlst_scheme_name(line[1])
            # incN -> IncN
            if re.search('incN', line[1], re.IGNORECASE):
                line[1] = re.sub('incN', 'IncN', line[1])
            # IncF: https://academic.oup.com/jac/article/65/12/2518/752763#83186263
            if re.search('incF', line[1], re.IGNORECASE):
                assert line[2] == '-', "IncF ST != -: {}({}): {}".format(line[1], line[2], ';'.join(line[3:]))
                # regular expressions
                re_k = r'\((?P<alleles>~?(\d+|-)\??(,~?(\d+|-)\??)*)\)' # to find alleles
                re_f = r'FI(?P<letter>(C|I|IS|IK|IY))\(' # to find FI(C/I/K/S/Y)
                re_a = r'FIA\(' # to find FIA
                re_b = r'FIB\(' # to find FIB
                # placeholders for alleles
                st_f = "-" # C/F/K/S/Y allele number
                st_l = "-" # C/F/K/S/Y letter if non-missing allele number
                st_a = "-" # A allele number
                st_b = "-" # B allele number
                mult_f = 0 # non-missing alleles for multiple letter from C/F/K/S/Y
                for locus in line[3:]:
                    # extract alleles
                    alleles = re.search(re_k, locus).group('alleles').split(',') # get and split alleles
                    # exact alleles: rm missing, with "~" or "?"
                    alleles_cleaned = [a for a in alleles if a != "-" and not re.search("~", a) and not re.search("\?", a)]
                    ge_zero = len(set(alleles_cleaned)) > 0 # at least one exact allele
                    mult_alleles = len(set(alleles_cleaned)) > 1 # multiple exact alleles
                    # which locus
                    if re.match(re_a, locus, re.IGNORECASE) and not mult_alleles:
                        if ge_zero:
                            st_a = alleles_cleaned[0]
                    elif re.match(re_b, locus, re.IGNORECASE) and not mult_alleles:
                        if ge_zero:
                            st_b = alleles_cleaned[0]
                    elif re.match(re_f, locus, re.IGNORECASE) and ge_zero:
                        mult_f += 1
                        # ambiguous: more than one loci with non-missing alleles
                        if mult_f > 1:
                            continue
                        # letter
                        st_l = re.match(re_f, locus).group('letter')
                        if st_l == "I":
                            st_l = "F"
                        else:
                            st_l = re.sub("I", "", st_l)
                        # allele
                        if mult_alleles:
                            continue
                        st_f = alleles_cleaned[0]
                if st_l == "-" and mult_f == 0: # no allele for C/F/K/S/Y
                    assert st_f == "-"
                    st_l = "F"
                elif mult_f > 1: # mult. loci with non-missing alleles
                    st_f = "-"
                    st_l = "-"
                line[2] = "{l}{f}:A{a}:B{b}".format(l=st_l, f=st_f, a=st_a, b=st_b)
            # save
            hits.append("{}({}):{}".format(line[1], line[2], ';'.join(line[3:])))
    return hits

##################################################
# MAIN
##################################################
ARGS = None
if __name__ == "__main__":
    # ARGS
    ARGS = get_arg_parser().parse_args()
    # Logger
    logger = My_logger(log_filename = ARGS.log, logger_name = "process_pmlst")
    logger = logger.get_logger()

    # PlasmidFinder
    pf = pandas.read_csv(ARGS.pf)
    pf = pf.loc[pf['sseqdb'] == 'plasmidfinder',]
    logger.info('Read in {} PlasmidFinder hits\n{}'.format(pf.shape[0], pf.head()))

    # pMLST
    pmlst_hits = []
    tmp_file = os.path.join('pmlst.tmp')
    tmp_fasta = os.path.join(tmp_file + '.fasta')
    logger.info('Run pMLST for each record with available scheme w.r.t. PlasmidFinder hits.')
    with open(ARGS.fasta, 'r') as ifile:
        for record in SeqIO.parse(ifile, 'fasta'):

            # found replicons
            record_replicons = list(pf.loc[pf['qseqid'] == record.id,'sseqid'])
            if len(record_replicons) == 0:
                continue
            else:
                SeqIO.write(record, open(tmp_fasta, 'w'), 'fasta') # create tmp FASTA

            used_schemes = set()
            for record_replicon in record_replicons:
                # set scheme
                scheme = None
                mapping = ARGS.map[0].replace('{', '').replace('}', '').replace('\'', '')
                mapping = [x for x in mapping.split(', ')]
                for m in mapping:
                    k = m.split(': ')[0]
                    v = m.split(': ')[1]
                    if re.match(k, record_replicon, re.IGNORECASE):
                        scheme = v
                        break
                if scheme is None or scheme in used_schemes: # no matching scheme/already done
                    continue
                else:
                    used_schemes.add(scheme)
                # run pmlst
                cmd = ARGS.cmd.format(scheme=scheme) + " {fasta} > {ofile}".format(fasta=tmp_fasta, ofile=tmp_file)
                run_cmd(cmd, split=False, shell=True)

                # process hits
                hits = process_pmlst_hits(f=tmp_file, 
                                          pmlst_db_path=os.path.join(ARGS.env_path, 'db/pubmlst'), 
                                          logger=logger)
                for i in range(0, len(hits)):
                    hits[i] = {
                        'ID': record.id,
                        'pmlst': hits[i]
                    }
                pmlst_hits += hits

    # add pMLST hits to PlasmidFinder hits
    pmlst_hits = pandas.DataFrame(pmlst_hits)

    # save
    pmlst_hits.to_csv(ARGS.ofile,index=False)

    # clean up
    os.remove(tmp_file)
    os.remove(tmp_fasta)