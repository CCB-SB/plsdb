# CODE TAKE FROM ABC-HUMI: https://github.com/gurevichlab/abc_humi/tree/main

import os
import sys
import csv
from collections import defaultdict
from pathlib import Path
from itertools import repeat
from tqdm import tqdm
import pandas as pd


BGC_classes = ['NRPS', 'Others', 'PKSI', 'PKS-NRP_Hybrids',
               'PKSother', 'RiPPs', 'Terpene']
MIBIG_PREF = 'BGC'
HUMI_PREF = 'HMBGC'


def get_representatives(gcf_trees_dir: Path, cutoff='0.30') -> dict:
    representatives = dict()
    for f in gcf_trees_dir.iterdir():
        if 'alignment' in f.name and f.name.startswith('GCF_c' + cutoff):
            fst_line = f.read_text().splitlines()[0]
            bgc_id = fst_line[1:].strip()
            # filename in format GCF_c0.30_02898_alignment.fasta
            family_id = str(int(f.name.split('_')[2]))
            representatives[family_id] = bgc_id

    return representatives


def get_clustering_data(class_dir: Path,
                        bgc_class: str,
                        cutoff_f='0.30',
                        cutoff_c='0.70') -> dict:
    clustering_file = class_dir / Path(f'{bgc_class}_clans_{cutoff_f}_{cutoff_c}.tsv')
    if clustering_file.exists():
        clustering_data = list(csv.DictReader(clustering_file.open('r'),
                                              delimiter='\t'))
    else:
        alt_clust_file = class_dir / Path(f'{bgc_class}_clustering_c{cutoff_f}.tsv')
        clustering_data = list(csv.DictReader(alt_clust_file.open('r'),
                                              delimiter='\t'))
        for row in clustering_data:
            row['Clan Number'] = row['Family Number']

    return clustering_data


def retrieve_big_scape_output_single_class(class_dir: Path,
                                           bgc_class: str,
                                           cutoff_f='0.30',
                                           cutoff_c='0.70') -> dict:
    representatives = get_representatives(class_dir / Path('GCF_trees'),
                                          cutoff=cutoff_f)
    clustering_data = get_clustering_data(class_dir, bgc_class,
                                          cutoff_f=cutoff_f,
                                          cutoff_c=cutoff_c)

    bgc_data = dict()
    for row in clustering_data:
        bgc_id = row['#BGC Name']
        family_id = row['Family Number']
        clan_id = row['Clan Number']

        # sometimes alignment file is absent but the family consists of several bgcs
        # perhaps that happens when there're only 2 bgcs in the family but I need to check
        # in this case I will carefully choose a representative later so that it is the same across all the bgc classes
        if family_id not in representatives:
            is_family_repr = False
            is_clan_repr = False
        else:
            is_family_repr = representatives[family_id] == bgc_id
            is_clan_repr = family_id == clan_id and is_family_repr

        bgc_data[bgc_id] = {'GCF_Type': bgc_class,
                            'GCC_ID': int(clan_id),
                            'GCF_ID': int(family_id),
                            'Is_representative_for_GCF': is_family_repr,
                            'Is_representative_for_GCC': is_clan_repr}
    return bgc_data


def retrieve_big_scape_output(bigscape_out_dir: Path,
                              cutoff_f='0.30', cutoff_c='0.70'):
    bgc_data = dict()
    for bgc_class in tqdm(BGC_classes):
        class_dir = bigscape_out_dir / Path(bgc_class)
        if not class_dir.exists():
            continue
        bgc_data[bgc_class] = retrieve_big_scape_output_single_class(class_dir,
                                                                     bgc_class,
                                                                     cutoff_f=cutoff_f,
                                                                     cutoff_c=cutoff_c)
    return bgc_data


def get_entries_in_families(big_scape_results, pref=MIBIG_PREF) -> dict:
    mibig_in_families = defaultdict(list)
    for bgc_class in big_scape_results:
        for bgc_id, bgc_data in big_scape_results[bgc_class].items():
            if bgc_id.startswith(pref):
                mibig_in_families[bgc_data['GCF_ID']].append(bgc_id)

    return mibig_in_families


def add_missing_representatives(big_scape_results):
    families = defaultdict(lambda: defaultdict(set))
    has_repr = set()
    for bgc_class in big_scape_results:
        for bgc_id, bgc_data in big_scape_results[bgc_class].items():
            families[bgc_data['GCF_ID']][bgc_data['GCF_Type']].add(bgc_id)
            if bgc_data['Is_representative_for_GCF']:
                has_repr.add(bgc_data['GCF_ID'])

    for family_id, family_data in families.items():
        if family_id not in has_repr:
            if any(len(family_data[bgc_class]) > 2
                   for bgc_class in family_data):
                print(f'Family {family_id} is big but has no representative!!!')
            shared_bgc = set.intersection(*(family_data[bgc_class]
                                            for bgc_class in family_data))
            if not shared_bgc:
                print('NO SHARED BGCs FOR', family_id)
                exit(0)
            rep = next(iter(shared_bgc))
            for bgc_class in big_scape_results:
                if family_id in big_scape_results[bgc_class]:
                    big_scape_results[bgc_class][rep]['Is_representative_for_GCF'] = True
                    if big_scape_results[bgc_class][rep]['GCC_ID'] == family_id:
                        big_scape_results[bgc_class][rep]['Is_representative_for_GCC'] = True

    return big_scape_results



    

if __name__ == '__main__':
    ##
    # MAIN
    ##

    # ARGS
    path_bigscape = os.path.join(snakemake.input.DIR, 'network_files')
    cutoff_gcc_gcf = snakemake.params.cutoff_gcf
    cutoff_gcc = snakemake.params.cutoff_gcc
    output_tsv = snakemake.output.tsv

    dirs = sorted(x for x in os.listdir(path_bigscape) if os.path.isdir(os.path.join(path_bigscape, x)))
    dir_latest = dirs[-1]
    path_latest = os.path.join(path_bigscape, dir_latest)

    bigscape_results = retrieve_big_scape_output(path_latest, cutoff_f=cutoff_gcc_gcf, cutoff_c=cutoff_gcc)
    bigscape_results = add_missing_representatives(bigscape_results)
    mibig_in_families = get_entries_in_families(bigscape_results, pref=MIBIG_PREF)
    humi_in_families = get_entries_in_families(bigscape_results, pref=HUMI_PREF)

    df = pd.DataFrame([
        {'BGC_ID': kk, 'GCF_Type': k, **vv, 
         'MIBiG': ','.join(mibig_in_families[vv['GCF_ID']]), 
         'ABCHuMi':  ','.join(humi_in_families[vv['GCF_ID']])}
        for k, v in bigscape_results.items()
        for kk, vv in v.items()
    ])

    representatives_gcf = df.loc[df['Is_representative_for_GCF'], ['GCF_ID', 'BGC_ID']].rename(columns={'BGC_ID':'Representative_GCF'})
    representatives_gcc = df.loc[df['Is_representative_for_GCC'], ['GCC_ID', 'BGC_ID']].rename(columns={'BGC_ID':'Representative_GCC'})
    
    df = pd.merge(df, representatives_gcf, how='outer', on='GCF_ID')
    df = pd.merge(df, representatives_gcc, how='outer', on='GCC_ID')

    df['GCF_ID'] = df['GCF_ID'].map(lambda x: f'GCF_{x:08d}')
    df['GCC_ID'] = df['GCC_ID'].map(lambda x: f'GCC_{x:07d}')
    
    mask = df['BGC_ID'].str.startswith(MIBIG_PREF) | df['BGC_ID'].str.startswith(HUMI_PREF)
    df = df.loc[~mask, ['BGC_ID', 'GCF_Type', 'GCC_ID', 'GCF_ID', 'MIBiG', 'ABCHuMi']]

    df.to_csv(output_tsv, sep='\t', index=False)