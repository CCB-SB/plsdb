#!/usr/bin/python

import os
import pickle
import pandas
import logging
import argparse
from tqdm import tqdm
from glob import glob
from Bio import Entrez
from Bio.SeqUtils import GC
from multiprocessing import Pool
import xml.etree.ElementTree as ET

from utils import *

##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    # Input/output
    parser.add_argument('--assemblies', '-a', help="Assembly tables", required=True, nargs="+")
    parser.add_argument('--sequences',  '-s', help="Sequence table", required=True)
    parser.add_argument('--emb',        '-e', help="Embedding", required=True)
    parser.add_argument('--locs',       '-l', help="Location information", required=False)
    parser.add_argument('--ofile',      '-o', help="Output file", required=True)
    parser.add_argument('--fa_dir', help="Directory containing genomic FASTA files", required=True)
    parser.add_argument('--cores', help="How many cores to use for parallel execution", type=int, default=1)
    parser.add_argument('--freshstart', help="Do not use checkpoints to restore data", action="store_true")
    parser.add_argument('--gmaps_api_keys', help="GoogleMaps API key", required=True, nargs="+")
    return parser

##################################################
# VARS
##################################################
# Relevant taxonomy ranks
TAX_RANKS = [
    'species',
    'genus',
    'family',
    'order',
    'class',
    'phylum',
    'superkingdom'
]

##################################################
# FUNCS
##################################################
def lindict2str(lin):
    """
    Lineage as list od dicts -> string
    """
    lin = ';'.join([
        '%s__%s__%s' % (
            re.sub('\s+', '_', l['Rank']),
            l['TaxId'],
            re.sub('\s+', '_', l['ScientificName'])
        ) for l in lin
    ])
    return lin

def proc_lindict(lin):
    """
    Lineage as list od dicts -> dict with relevant info
    """
    # primary taxon
    new_lin = {
        'taxon_rank': lin['Rank'],
        'taxon_id':   lin['TaxId'],
        'taxon_name': lin['ScientificName']
    }
    # lineage + primary taxon
    lin = lin['LineageEx'] + [{
        'Rank': lin['Rank'],
        'TaxId': lin['TaxId'],
        'ScientificName': lin['ScientificName']
    }]
    # lineage as str
    new_lin['lineage'] =  lindict2str(lin)
    # relevant ranks
    for rank in TAX_RANKS:
        rank_taxon = [t for t in lin if t['Rank'] == rank]
        assert len(rank_taxon) <= 1
        if len(rank_taxon) == 1:
            new_lin['taxon_%s_id' % rank]   = rank_taxon[0]['TaxId']
            new_lin['taxon_%s_name' % rank] = rank_taxon[0]['ScientificName']
        else:
            new_lin['taxon_%s_id' % rank]   = None
            new_lin['taxon_%s_name' % rank] = None
    return new_lin

def prep_ids(ids):
    """
    Prepare list of IDs for E-utils
    """
    ids = sorted(list(set(ids)))
    return ids, ','.join([str(i) for i in ids])

def get_lineages(ids):
    """
    Extract taxonomic lineages for given taxonomy IDs
    returns a dict with (<tax. ID>: <dict>)
    """
    ids, ids_str = prep_ids(ids)
    lineages = {}
    for i, record in enumerate(Entrez.parse(Entrez.efetch(db="taxonomy", id=ids_str))):
        if str(ids[i]) != record['TaxId']:
            logging.info('WARNING: Match for TaxId %s is %s' % (str(ids[i]), record['TaxId']))
        lineages[ids[i]] = proc_lindict(record)
    return lineages

def get_seqsmeta(ids):
    """
    Extract sequence information from nuccore
    returns a dict with (<seq. ID>: <dict>)
    """
    ids, ids_str = prep_ids(ids)
    seqs_meta = {}
    records   = [r for r in Entrez.parse(Entrez.efetch(db="nuccore", id=ids_str, retmode="xml"))]
    for i, record in enumerate(records):
        assert ids[i] == record['GBSeq_accession-version'], 'Not matching IDs %s and %s' % (ids[i], record['GBSeq_accession-version'])
        seqs_meta[ids[i]] = {
            'strandedness': record['GBSeq_strandedness'],
            'moltype': record['GBSeq_moltype'],
            'topology': record['GBSeq_topology'],
            'create_date': record['GBSeq_create-date'],
            'update_date': record['GBSeq_update-date']
        }
    return seqs_meta

def get_biosamplemeta(ids):
    """
    Extract BioSample information
    Need XML parser as Entrez.read/parse does not work for BioSample
    returns a dict with (<BioSample ID>: <dict>)
    """
    ids, ids_str = prep_ids(ids)
    handle = Entrez.efetch(db="biosample", id=ids_str, retmode="xml")
    records = ET.parse(handle).getroot()
    bs_meta = {}
    keys = ['geo_loc_name', 'lat_lon', 'isolation_source', 'specific_host', 'host_disease', 'host_tissue_sampled', 'host_description', 'collection_date']
    for record in records:
        # BioSample ID
        assert record.tag == "BioSample", record.tag
        record_id = record.attrib['accession']
        bs_meta[record_id] = {}
        # Get data
        for c in record:
            if c.tag == 'Attributes':
                for cc in c:
                    if cc.attrib['attribute_name'] in keys:
                        bs_meta[record_id][cc.attrib['attribute_name']] = cc.text
    return bs_meta

def get_seqsstats(ffile):
    stats = dict()
    for record in fasta_records(ffile):
        stats[record.id] = {
            'sequence_length': len(record.seq),
            'sequence_gc': GC(str(record.seq))
        }
    return stats

def collect_lineages(pool, ids, chunk_num):
    logging.info('Retrieve taxonomy lineages')
    tmp = pool.map(
        get_lineages,
        [chunk for chunk in chunk_list(set(ids), chunk_num)]
    )
    # lists to dict
    lineages = {}
    for d in tmp:
        lineages.update(d)
    # dict to dataframe
    lineages = pandas.DataFrame.from_dict(lineages, orient="index")
    # lineages = lineages.assign(query_id=lineages.index)
    return lineages

def collect_seqsmeta(pool, ids, chunk_num):
    logging.info('Retrieve sequence meta data')
    tmp = pool.map(
        get_seqsmeta,
        [chunk for chunk in chunk_list(set(ids), chunk_num)]
    )
    # lists to dict
    seq_meta = {}
    for d in tmp:
        seq_meta.update(d)
    # dict to dataframe
    seq_meta = pandas.DataFrame.from_dict(seq_meta, orient="index")
    return seq_meta

def collect_biosamplemeta(pool, ids, chunk_num):
    logging.info('Retrieve BioSample meta data')
    tmp = pool.map(
        get_biosamplemeta,
        [chunk for chunk in chunk_list(set(ids), chunk_num)]
    )
    # lists to dict
    bs_meta = {}
    for d in tmp:
        bs_meta.update(d)
    # dict to dataframe
    bs_meta = pandas.DataFrame.from_dict(bs_meta, orient="index")
    return bs_meta

def collect_seqsstats(pool, idir, df):
    # find files
    fa_files = sorted(glob("%s/*_genomic.fna.gz" % idir))
    logging.info("Found %d genomic FASTA files" % len(fa_files))
    # collect data
    tmp = pool.map(get_seqsstats, fa_files)
    # merge results
    stats = {}
    for d in tmp:
        stats.update(d)
    # add to df
    df = df.assign(sequence_length=[stats[k]["sequence_length"] for k in df.index])
    df = df.assign(sequence_gc=[stats[k]["sequence_gc"] for k in df.index])
    return df

##################################################
# MAIN
##################################################

if __name__ == "__main__":
    # Logger setup
    setup_logger()

    # Args
    args = get_arg_parser().parse_args()

    # Checkpoints
    cp = {
        "lins": "%s_lins.pck" % os.path.splitext(args.ofile)[0],
        "meta": "%s_meta.pck" % os.path.splitext(args.ofile)[0],
        "bios": "%s_bios.pck" % os.path.splitext(args.ofile)[0],
        "locs": "%s_locs.pck" % os.path.splitext(args.ofile)[0],
        "seqs-fa": "%s_seqs-fa.pck" % os.path.splitext(args.ofile)[0],
    }

    # E-utils settings
    Entrez.email = "Your.Name.Here@example.org"

    # Open pool
    pool = Pool(args.cores)

    ##############################################
    # Sequences
    seqs = pandas.read_csv(filepath_or_buffer=args.sequences, sep='\t', header=0)
    seqs.set_index(keys='sequence_accession', drop=False, inplace=True, verify_integrity=True)
    seqs.drop(labels='asm_name', axis='columns', inplace=True) # NOTE derived from FASTA file name, not needed here
    logging.info("Sequences: %d entries" % seqs.shape[0])

    ##############################################
    # Assemblies
    asms = read_ncbi_tables(args.assemblies)
    asms.set_index(keys='assembly_accession', drop=False, inplace=True, verify_integrity=True)
    logging.info("Assemblies: %d entries" % asms.shape[0])

    ##############################################
    # Embedding
    emb = pandas.read_csv(
        filepath_or_buffer=args.emb,
        sep='\t',
        header=0,
        index_col=0
    )

    ##############################################
    # Taxonomic lineages
    if args.freshstart or not os.path.exists(cp["lins"]):
        lineages = collect_lineages(pool=pool, ids=asms['taxid'], chunk_num=args.cores*10)
        with open(cp["lins"], "wb") as ofile:
            pickle.dump(lineages, ofile)
    else:
        logging.info('Restoring checkpoint data from %s' % cp["lins"])
        with open(cp["lins"], "rb") as ifile:
            lineages = pickle.load(ifile)
    logging.info("Lineages: %d entries" % lineages.shape[0])

    ##############################################
    # Meta data
    if args.freshstart or not os.path.exists(cp["meta"]):
        seq_meta = collect_seqsmeta(pool=pool, ids=seqs['sequence_accession'], chunk_num=args.cores*100)
        with open(cp["meta"], "wb") as ofile:
            pickle.dump(seq_meta, ofile)
    else:
        logging.info('Restoring checkpoint data from %s' % cp["meta"])
        with open(cp["meta"], "rb") as ifile:
            seq_meta = pickle.load(ifile)
    logging.info("Meta data: %d entries" % seq_meta.shape[0])

    ##############################################
    # Sequence statistics
    if args.freshstart or not os.path.exists(cp["seqs-fa"]):
        seq_meta = collect_seqsstats(pool=pool, idir=args.fa_dir, df=seq_meta)
        with open(cp["seqs-fa"], "wb") as ofile:
            pickle.dump(seq_meta, ofile)
    else:
        logging.info('Restoring checkpoint data from %s' % cp["seqs-fa"])
        with open(cp["seqs-fa"], "rb") as ifile:
            seq_meta = pickle.load(ifile)
    logging.info('Added sequence-derived statistics')

    ##############################################
    # BioSample meta data
    if args.freshstart or not os.path.exists(cp["bios"]):
        bs_meta = collect_biosamplemeta(pool=pool, ids=set(asms['biosample']), chunk_num=args.cores*100)
        with open(cp["bios"], "wb") as ofile:
            pickle.dump(bs_meta, ofile)
    else:
        logging.info('Restoring checkpoint data from %s' % cp["bios"])
        with open(cp["bios"], "rb") as ifile:
            bs_meta = pickle.load(ifile)
    logging.info("BioSample meta data: %d entries" % bs_meta.shape[0])

    ##############################################
    # Parse BioSample location data
    bs_meta['loc_lat'] = None
    bs_meta['loc_lng'] = None
    if args.freshstart or not os.path.exists(cp["locs"]):
        if args.locs is not None:
            logging.info('Reading location mapping from %s' % args.locs)
            with open(args.locs, "rb") as ifile:
                locs = pickle.load(ifile)
        else:
            locs = {}
    else:
        logging.info('Restoring checkpoint data from %s' % cp["locs"])
        with open(cp["locs"], "rb") as ifile:
            locs = pickle.load(ifile)
    current_key = args.gmaps_api_keys.pop()
    for i in bs_meta.index:
        # location name and coordinates
        l_n = bs_meta.loc[i,'geo_loc_name']
        l_c = bs_meta.loc[i,'lat_lon']
        # avoid nulls
        if pandas.isnull(l_n) or l_n is None:
            l_n = ''
        if pandas.isnull(l_c) or l_c is None:
            l_c = ''
        # location is not already saved
        if (l_n, l_c) not in locs:
            logging.info('Retrieving coordinates for location: \"%s\", \"%s\"' % (l_n, l_c))
            try_new_key = True
            set_new_key = False
            while try_new_key:
                try:
                    loc = parse_location(loc_name=l_n, loc_coords=l_c, api_key=current_key)
                    # save location
                    locs[(l_n,l_c)] = loc
                    with open(cp["locs"], "wb") as ofile:
                        pickle.dump(locs, ofile)
                    # new key worked
                    try_new_key = False
                except Exception as e:
                    if set_new_key:
                        logging.error('Error while retrieving coordinates for location with new key: %s, %s' % (l_n, l_c))
                        raise(e)
                    else:
                        logging.error('Error while retrieving coordinates for location: %s, %s' % (l_n, l_c))
                        assert len(args.gmaps_api_keys) > 0, 'No new API keys'
                        current_key = args.gmaps_api_keys.pop()
                        try_new_key = True
                        set_new_key = True
        # if (l_n, l_c) not in locs or locs[(l_n,l_c)] is None: # TODO: rm
        if locs[(l_n,l_c)] is None:
            bs_meta.loc[i,'loc_lat'] = None
            bs_meta.loc[i,'loc_lng'] = None
        else:
            bs_meta.loc[i,'loc_lat'] = locs[(l_n,l_c)]['lat']
            bs_meta.loc[i,'loc_lng'] = locs[(l_n,l_c)]['lng']
    if args.freshstart or not os.path.exists(cp["locs"]):
        with open(cp["locs"], "wb") as ofile:
            pickle.dump(locs, ofile)

    # Close pool
    pool.close(); pool.join()

    ##############################################
    # Process/merge
    # Filter assemblies - keep only the relevant ones
    asms = asms[asms['assembly_accession'].isin(seqs['assembly_accession'])]
    logging.info("Assemblies: %d entries were kept" % asms.shape[0])

    # Add taxonomy to assemblies
    asms = pandas.merge(
        left=asms,
        right=lineages,
        how='left',
        left_on='taxid',
        right_index=True,
        sort=False
    )
    assert len([i for i in asms['taxon_id'] if i is None]) == 0, 'Not all assemblies have a taxonomy match'
    logging.info("Assemblies: taxonomy info was added - %d entries" % asms.shape[0])

    # Add BioSample meta data to assemblies
    asms = pandas.merge(
        left=asms,
        right=bs_meta,
        how='left',
        left_on='biosample',
        right_index=True,
        sort=False
    )
    logging.info("Assemblies: BioSample meta data was added - %d entries" % bs_meta.shape[0])

    # Merge seqs and seqs meta
    seqs = pandas.merge(
        left=seqs,
        right=seq_meta,
        how='left',
        left_index=True,
        right_index=True,
        sort=False
    )
    assert len([i for i in seqs['strandedness'] if i is None]) == 0, 'Not all sequences have a meta data match'
    logging.info("Sequences: meta data was added - %d entries" % seqs.shape[0])

    # Add embedding
    seqs = pandas.merge(
        left=seqs,
        right=emb,
        how='left',
        left_index=True,
        right_index=True,
        sort=False
    )
    assert len([i for i in seqs['D1'] if i is None]) == 0, 'Not all sequences have embedding coordinates'
    logging.info("Sequences: embedding was added - %d entries" % seqs.shape[0])

    ##############################################
    # Master table
    master = pandas.merge(
        left=seqs,
        right=asms,
        how='outer',
        left_on="assembly_accession",
        right_index=True,
        sort=False
    )
    assert len([i for i in master['asm_name'] if i is None]) == 0, 'Not all sequences have an assembly match'
    logging.info("Master table was created - %d entries" % master.shape[0])
    # drop some columns
    drop_cols = ["assembly_accession_x", "assembly_accession_y"]
    master.drop(labels=[h for h in drop_cols if h in master.columns], axis='columns', inplace=True)

    # Save
    master.to_csv(path_or_buf=args.ofile, sep='\t', header=True, index=False, index_label=False)
