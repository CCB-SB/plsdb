    
#! /usr/bin/env/python
# coding: utf-8

## -----
# Retrival of NCBI data
# Author: Unknown person and modified by Alejandra Gonzalez
# Last modified:
## -----

###########
# IMPORTS #
###########
import argparse
from utils_my_logger import My_logger
from utils_decorators import timer, debuger

import logging
logg = logging.getLogger("pmlst_download")
###########
# FUNCTIONS #
###########

def get_locus_info(locus_path, scheme_dir, scheme_name, logger=None):
    import requests
    from os.path import join

    locus_info = requests.get(locus_path)
    locus_info = locus_info.json()
    locus = locus_info['id']
    if locus_info['alleles_fasta']:
        seqs = requests.get(locus_info['alleles_fasta'])
        out = join(scheme_dir, f'{locus}.tfa')
        with open(out, 'w+') as ofile:
            ofile.write(seqs.text)
    else:
        logger.warning(f'Scheme {scheme_name}, locus {locus}: no alleles')

@timer(my_logger=logg)
@debuger(my_logger=logg)
def download_pmlst_scheme_profiles(scheme_name, scheme_url, scheme_dir, scheme_profiles, logger = None):
    """
    Download pMLST scheme profiles for given scheme
    Will also check the formatting s.t. the files can be used by the "mlst" tool
    :param scheme name: name of the scheme
    :param scheme_url: URL to be used
    :param scheme_dir: where the allele FASTA files (*.tfa) are stored
    :param scheme_profiles: output file name for the profiles
    """
    import pandas as pd
    import os
    import re
    import requests
    from glob import glob

    profiles = requests.get(scheme_url + '/profiles_csv')
    
    # save profiles
    if profiles.status_code == 404: # no profiles
        allele_files = sorted(glob('%s/*.tfa' % scheme_dir))
        loci = [ re.sub('\.tfa$', '', os.path.basename(tfa)) for tfa in allele_files ]
        # profiles file with a dummy entry
        with open(scheme_profiles, 'w') as ofile:
            # header
            ofile.write('ST\t{}\n'.format('\t'.join(loci)))
            # dummy entry: ST + loci
            ofile.write('\t'.join( ['1']*( len(loci)+1 ) ))
        # empty file saying that the profiles file contains dummy data
        with open(scheme_profiles + '.dummy', 'w') as ofile:
            ofile.write('dummy')
    else: # save profiles
        with open(scheme_profiles, 'w') as ofile:
            ofile.write(re.sub('\t\n', '\n', profiles.text))
    assert os.path.exists(scheme_profiles), 'No profiles file for scheme {}'.format(scheme_name)

    # check formatting
    df = pd.read_csv(scheme_profiles, sep='\t', header=0)
    assert df.shape[0] > 0, 'Empty profiles file for scheme {}'.format(scheme_name)

    # 1st column must be "ST"
    if df.columns[0] != "ST":
        df.columns = ['ST'] + list(df.columns)[1:]

    # all STs must be integers
    if any([not re.fullmatch(r'\d+', str(st)) for st in list(df['ST'])]):
        logger.warning('Scheme {}: not all ST values are integers'.format(scheme_name))
        # copy of original values
        df['oldST'] = df['ST'].copy()
        # number STs from 1 to N
        df['ST'] = list(range(1, df.shape[0]+1))
        # save old values
        df.to_csv(scheme_profiles + '.old', sep='\t', index=False, index_label=False)

    # save
    cols = [col for col in list(df.columns) if col != 'oldST']
    df[cols].to_csv(scheme_profiles, sep='\t', index=False, index_label=False)
    return

@timer(my_logger=logg)
@debuger(my_logger=logg)
def download_pmlst_scheme_alleles(scheme_name, scheme_url, scheme_dir, logger = None, cores=1):
    """
    Download pMLST scheme allele sequences
    :param scheme_name: name of the scheme
    :param: scheme_url: URL to be used
    :param scheme_dir: output directory
    """
    from os.path import join
    import requests
    import multiprocessing as mp

    # get loci
    loci = requests.get(scheme_url.replace('isolates', 'seqdef')).json()['loci']
    
    # download alleles for each locus
    args = [(locus_path, scheme_dir, scheme_name, logger) for locus_path in loci]
    
    with mp.Pool(processes = cores) as pool:
      try:
          pool.starmap(get_locus_info, args)
      except Exception as e:
          logger.error("Something when wrong Downloading allelese for each locus")
          pool.terminate()
          raise e

    return

def proc_mlst_scheme_name(scheme_name):
    """
    Process scheme name, e.g. to be used as directory name
    """
    import re
    return re.sub('/', '_', re.sub('\s+', '__', scheme_name))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Obtain plasmid NCBI query ",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Mandatory argument
    parser.add_argument("--db", required = True, 
                        help="Terms to include in the esearch query")
    parser.add_argument("--url", required = True, 
                        help="NCBI API KEY")
    parser.add_argument("--url-schemes", required = True,
                        help="Valid email for NCBI query" )
    parser.add_argument("--outdir", default = "results", required = False,
                        help="Output directory" )
    parser.add_argument("--outfile", default = "results.fna", required = False,
                        help="Output file" )
    parser.add_argument("--log", default = "log.log", required = False)
    parser.add_argument("--cores", default = 1, type=int, required = False)

    args = parser.parse_args()
  
    # reference: https://github.com/kjolley/BIGSdb.git: scripts/rest_examples/python/download_alleles.py
    import requests
    from os.path import join
    import os
    import shutil
    
    
    # Logging
    logger = My_logger(log_filename = args.log, logger_name = "pmlst_download")
    logger = logger.get_logger()

    # get schemes
    schemes = requests.get(args.url_schemes)
    assert schemes.status_code != 404, f'Invalid URL {args.url_schemes}'
    schemes = schemes.json()['schemes']
    logger.info(f'There are {len(schemes)} pMLST schemes')

    # get loci
    for scheme in schemes:
        # scheme name and URL
        scheme_name = scheme['description']
        scheme_url  = scheme['scheme']
        scheme_name2 = proc_mlst_scheme_name(scheme_name)
        logger.info(f'Scheme {scheme_name} (directory: {scheme_name2}): {scheme_url}')

        # Create specific directory for each scheme
        scheme_dir = join(args.outdir, scheme_name2)
        if (os.path.isdir(scheme_dir)):
            shutil.rmtree(scheme_dir)
        
        os.mkdir(scheme_dir)

        # Scheme profile file
        scheme_profiles = join(scheme_dir, f'{scheme_name2}.txt')

        # alleles
        download_pmlst_scheme_alleles(scheme_name=scheme_name, scheme_url = scheme_url,
                                      scheme_dir = scheme_dir, logger = logger, cores=args.cores)

        # profiles
        profiles_url = f"{args.url}/db/{args.db}/schemes/{scheme_url.split('/')[-1]}"
        download_pmlst_scheme_profiles(scheme_name = scheme_name, scheme_url = profiles_url, 
                                       scheme_dir = scheme_dir, scheme_profiles = scheme_profiles,
                                       logger = logger)

    with open(args.outfile, 'w') as ofile:
        ofile.write('done')