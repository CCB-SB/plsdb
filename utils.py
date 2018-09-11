#!/usr/bin/python

import re
import logging
import pandas

##################################################
# MISC
##################################################
def setup_logger(log_level=logging.DEBUG):
    """Logger"""
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%y:%m:%d %H:%M:%S'
    )

def split_list(values, size):
    """
    Split given list into chunks of given size
    """
    values = list(values)
    for i in range(0, len(values), size):
        yield values[i:(i + size)]

def run_cmd(cmd):
    """
    Run given CMD
    :param cmd: CMD (str)
    :return: cmd (str), status (int), stdout + stderr (str)
    """
    import subprocess
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    p_stdout = p.stdout.read().decode()
    p_comm   = p.communicate()[0]
    p_status = p.returncode
    return cmd, p_status, p_stdout

def mkdir(dpath, p=True):
    """
    Mimick "mkdir -p"
    """
    import os
    if p:
        if not os.path.exists(dpath):
            os.makedirs(dpath)
    else:
        os.makedirs(dpath)

def str2timestamp(ts, ts_format='%Y-%m-%d %H:%M:%S'):
    """
    Time stamp string to time stamp object
    """
    import datetime
    return datetime.datetime.strptime(ts, ts_format)

##################################################
# TOOLS
##################################################
def run_epost_split(df_file, ofile, header, cmd, df_col, split_size, split_str, **kwargs):
    """
    To run epost queries by splitting the input IDs into chunks (as the max. number of IDs is limited).
    Checks for unexpected IDs
    :param df_file: Data table containing needed IDs
    :param ofile: Output file
    :param header: Output file header
    :param cmd: epost CMD (con contain subsequent calls of other tools, e.g. efilter, xtract etc.)
    :param df_col: column name containing the IDs in the given table
    :param split_size: How large the ID chunks should be
    :param split_str: If ID-column values should be split by given string to get IDs, i.e. a row could containg a list of IDs
    :param **kwargs: key-word parameters for the CMD
    """
    import pandas, tqdm

    logging.info('CMD: {}'.format(cmd))

    ids = pandas.read_csv(df_file, sep='\t', header=0, dtype=str)[df_col]
    logging.info('There are {} IDs'.format(len(ids)))
    ids.fillna(value='', inplace=True)
    ids = ids.loc[ids != ""]
    ids = set(ids.values)
    if split_str is not None and split_str != "":
        ids_ = set()
        for id_ in ids:
            ids_ = ids_.union(id_.split(split_str))
        ids = ids_
    logging.info('There are {} unique IDs'.format(len(ids)))

    with open(ofile, 'w') as ohandle:
        ohandle.write(header + '\n')
        for chunk in tqdm.tqdm(split_list(ids, split_size)):
            cmd_c = cmd.format(**kwargs, ids=','.join(chunk))
            cmd_c, cmd_s, cmd_o = run_cmd(cmd_c)
            assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd_c, cmd_s, cmd_o)
            # checks
            if re.search('elink', cmd):
                for line in cmd_o.split('\n'):
                    if line == '':
                        continue
                    if re.search('\t', line):
                        i_n, i_s = line.split('\t')
                    else:
                        i_n = line.strip()
                        i_s = None
                    assert i_n in ids, 'Unexpected nuccore ID {}'.format(i_n)
                    # assert i_s is None or (not re.search(';', i_s)), 'Multiple hits for {}: {}'.format(i_n, i_s)
            else:
                found = [l.split('\t')[0] for l in cmd_o.split('\n') if l != '']
                not_found = list(set(chunk).difference(found)) # not found IDs
                une_found = list(set(found).difference(chunk)) # unexpected IDs, i.e. should not be searched
                # assert len(not_found) == 0, 'Not found: {}: {}, ...'.format(len(not_found), ', '.join([str(s) for s in not_found[0:5]]))
                assert len(une_found) == 0, 'Not expected: {}: {}, ...'.format(len(une_found), ', '.join([str(s) for s in une_found[0:5]]))
                if len(not_found) > 0:
                    logging.warn('Could not find: %s' % ', '.join([str(s) for s in not_found]))
            # save
            ohandle.write(cmd_o)

def run_blastn_check(acc, obname, main_fasta, blastn_cmd, blastn_bin, blastn_header, blastn_pident):
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
    acc_fasta = '%s.%s.fna' % (obname, acc)
    acc_ofile = '%s.%s.tsv' % (obname, acc)
    with open(main_fasta, 'r') as ifile:
        for record in SeqIO.parse(ifile, 'fasta'):
            if record.id == acc:
                SeqIO.write(record, acc_fasta, 'fasta')
                break
    assert os.path.exists(acc_fasta), 'No file {}'.format(acc_fasta)
    cmd = blastn_cmd.format(
        bin=blastn_bin,
        fasta=acc_fasta,
        output=acc_ofile,
        header=blastn_header,
        pident=blastn_pident
    )
    logging.info('START {}'.format(acc))
    cmd, cmd_s, cmd_o = run_cmd(cmd)
    assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd, cmd_s, cmd_o)
    logging.info('END {}'.format(acc))
    return

def proc_mlst_scheme_name(scheme_name):
    import re
    return re.sub('/', '_', re.sub('\s+', '__', scheme_name))

def reproc_mlst_scheme_name(scheme_name):
    import re
    return re.sub('_', '/', re.sub('__', ' ', scheme_name))

def download_pmlst_scheme_alleles(scheme_name, scheme_url, scheme_dir):
    """
    TODO
    """
    import os
    import requests

    # get loci
    loci = requests.get(scheme_url.replace('isolates', 'seqdef')).json()['loci']
    # download alleles for each locus
    for locus_path in loci:
        locus_info = requests.get(locus_path)
        locus_info = locus_info.json()
        locus = locus_info['id']
        if locus_info['alleles_fasta']:
            seqs = requests.get(locus_info['alleles_fasta'])
            with open(os.path.join(scheme_dir, '%s.tfa' % locus), 'w') as ofile:
                ofile.write(seqs.text)
        else:
            logging.warning('Scheme {}, locus {}: no alleles'.format(scheme_name, locus))
    return

def download_pmlst_scheme_profiles(scheme_name, scheme_url, scheme_dir, scheme_profiles):
    """
    TODO:
    """
    import os
    import re
    import requests
    from glob import glob

    # r = requests.get(params.url +  '/db/' + params.db + '/schemes/' + str(scheme_id) + '/profiles_csv')
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
    df = pandas.read_csv(scheme_profiles, sep='\t', header=0)
    assert df.shape[0] > 0, 'Empty profiles file for scheme {}'.format(scheme_name)

    # 1st column must be "ST"
    if df.columns[0] != "ST":
        df.columns = ['ST'] + list(df.columns)[1:]

    # all STs must be integers
    if any([not re.fullmatch(r'\d+', str(st)) for st in list(df['ST'])]):
        logging.warning('Scheme {}: not all ST values are integers'.format(scheme_name))
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

##################################################
# LOCATIONS
##################################################
# Parsing location information from NCBI (BioSample)

def load_locs(loc_file):
    import pandas, os
    locs = None
    if os.path.exists(loc_file):
        locs = pandas.read_csv(loc_file, sep='\t', header=0)
        locs.set_index('location', drop=False, inplace=True, verify_integrity=True)
        logging.info('Loaded locations from {}\n{}'.format(loc_file, locs.head()))
    else:
        logging.info('No known locations to load.')
    return locs

def update_locs(locs, dict):
    import pandas
    if locs is None:
        locs = pandas.DataFrame([dict])
    else:
        locs = locs.append(pandas.DataFrame([dict]))
    locs.set_index('location', drop=False, inplace=True, verify_integrity=True)
    return locs

def save_locs(locs, ofile):
    assert locs is not None
    locs.to_csv(ofile, sep='\t', header=True, index=False, index_label=False)
    return

# Strings used for missing locations (to avoid unneccessary queries)
location_missing = [
    "",
    "-",
    "na", "n/a", "n.a.",
    "missing", "none",
    "not applicable", "not available", "not collected", "not determined", "not recorded",
    "unavailable", "unknown", "unspecified",
]

def loc_is_missing(loc_str):
    for na_str in location_missing:
        if re.fullmatch(na_str, loc_str, flags=re.IGNORECASE):
            return True
    return False

def handle_loc_exceptions(loc_str):
    import re
    if re.search(r'BIFSCo\s+Region\s+\d?', loc_str, flags=re.IGNORECASE) and re.search(r'US|USA|United States', loc_str, flags=re.IGNORECASE):
        loc_str = 'USA'
    elif re.search(r'Australia', loc_str, flags=re.IGNORECASE) and re.search(r'sealake', loc_str, flags=re.IGNORECASE):
        loc_str = re.sub('sealake', 'Sea Lake', loc_str)
    elif re.search(r'South Korea', loc_str, flags=re.IGNORECASE) and re.search(r'Yeo-?su.+sediment', loc_str, flags=re.IGNORECASE):
        loc_str = 'South Korea, Yeosu'
    elif re.search(r'South Korea', loc_str, flags=re.IGNORECASE) and re.search(r'Shnan-gun', loc_str, flags=re.IGNORECASE):
        loc_str = 'South Korea'
    elif re.search(r'South Korea', loc_str, flags=re.IGNORECASE) and re.search(r'Geoje\s+Island', loc_str, flags=re.IGNORECASE):
        loc_str = 'South Korea,Geoje'
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Harbin\s+Veterinary\s+Research\s+Institute', loc_str, flags=re.IGNORECASE):
        loc_str = 'China,Harbin'
    elif re.search(r'Germany', loc_str, flags=re.IGNORECASE) and re.search(r'Black\s+Forest', loc_str, flags=re.IGNORECASE):
        loc_str = 'Germany,Black Forest'
    elif re.search(r'Sweden', loc_str, flags=re.IGNORECASE) and re.search(r'Kosterfjord', loc_str, flags=re.IGNORECASE):
        loc_str = 'Sweden,Kosterfjorden'
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Bo\s*Hai.+Panjin', loc_str, flags=re.IGNORECASE):
        loc_str = 'China,Panjin'
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Sunitezuoqi', loc_str, flags=re.IGNORECASE):
        loc_str = 'China' # location unknown even in Google Maps
    elif re.search(r'Brazil', loc_str, flags=re.IGNORECASE) and re.search(r'Rondonia', loc_str, flags=re.IGNORECASE):
        loc_str = 'Brazil,Rondonia'
    elif re.search(r'South\s+Korea', loc_str, flags=re.IGNORECASE) and re.search(r'Tae-an', loc_str, flags=re.IGNORECASE):
        loc_str = 'South Korea,Taean'
    elif re.search(r'South\s+Africa', loc_str, flags=re.IGNORECASE) and re.search(r'Kruger\s+National\s+Park', loc_str, flags=re.IGNORECASE):
        loc_str = 'South Africa,Kruger National Park'
    elif re.search(r'Japan', loc_str, flags=re.IGNORECASE) and re.search(r'Himeji.*Univrersity\s+of\s+Hyogo', loc_str, flags=re.IGNORECASE):
        loc_str = 'Japan,Himeji'
    elif re.search(r'USA', loc_str, flags=re.IGNORECASE) and re.search(r'Baytown.*Burnet\s+Shores', loc_str, flags=re.IGNORECASE):
        loc_str = 'USA,Baytown'
    elif re.search(r'Thailand', loc_str, flags=re.IGNORECASE) and re.search(r'Nakhornrachisma', loc_str, flags=re.IGNORECASE):
        loc_str = re.sub('Nakhornrachisma', 'Nakhon Ratchasima', loc_str, flags=re.IGNORECASE)
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Xinjiang\s+Province.*desert', loc_str, flags=re.IGNORECASE):
        loc_str = 'China,Xinjiang'
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Eastern\s+Hubei\s+Province', loc_str, flags=re.IGNORECASE):
        loc_str = 'China,Hubei'
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Chenmai\s+qiaotou', loc_str, flags=re.IGNORECASE):
        loc_str = 'China,Qiaotouzhen'
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Hongyuan\s+Prairie', loc_str, flags=re.IGNORECASE):
        loc_str = 'China,Hongyuan'
    elif re.search(r'Brazil', loc_str, flags=re.IGNORECASE) and re.search(r'Tupasi', loc_str, flags=re.IGNORECASE):
        loc_str = 'Brazil'
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Xinjiang\s+Province', loc_str, flags=re.IGNORECASE):
        loc_str = 'China,Xinjiang'
    elif re.search(r'USA', loc_str, flags=re.IGNORECASE) and re.search(r'Bear\s+River\s+Refuge', loc_str, flags=re.IGNORECASE):
        loc_str = 'USA,Bear River Migratory Bird Refuge' # the only matching refuge location
    elif re.search(r'USA', loc_str, flags=re.IGNORECASE) and re.search(r'Ogden\s+Bay\s+Refuge', loc_str, flags=re.IGNORECASE):
        loc_str = 'USA,Ogden'
    elif re.search(r'Vietnam', loc_str, flags=re.IGNORECASE) and re.search(r'Do Xongpha', loc_str, flags=re.IGNORECASE):
        loc_str = 'Vietnam' # could not find the associated location
    elif re.search(r'India', loc_str, flags=re.IGNORECASE) and re.search(r'Haiderabad/Hindustan', loc_str, flags=re.IGNORECASE):
        loc_str = re.sub('Haiderabad/Hindustan', 'Haiderabad', loc_str, flags=re.IGNORECASE)
    elif re.search(r'Democratic\s+Republic\s+of\s+the\s+Congo', loc_str, flags=re.IGNORECASE) and re.search(r'Iturie\s+province', loc_str, flags=re.IGNORECASE):
        loc_str = 'Democratic Republic of the Congo,Ituri'
    elif re.search(r'Germany', loc_str, flags=re.IGNORECASE) and re.search(r'Bruschal', loc_str, flags=re.IGNORECASE):
        loc_str = 'Germany,Bruchsal'
    elif re.search(r'Antarctica', loc_str, flags=re.IGNORECASE) and re.search(r'Mount\s+Rittmann', loc_str, flags=re.IGNORECASE):
        loc_str = 'Antarctica'
    elif re.search(r'Canada', loc_str, flags=re.IGNORECASE) and re.search(r'University\s+of\s+British\s+Columbia', loc_str, flags=re.IGNORECASE):
        loc_str = 'Canada,Vancouver'
    elif re.search(r'Brazil', loc_str, flags=re.IGNORECASE) and re.search(r'Ribeirao\s+Preto.*Sao\s+Paulo\s+State', loc_str, flags=re.IGNORECASE):
        loc_str = 'Brazil,Ribeirao Preto'
    elif re.search(r'Japan', loc_str, flags=re.IGNORECASE) and re.search(r'Niigata.*Nagakura', loc_str, flags=re.IGNORECASE):
        loc_str = 'Japan,Nagakura'
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Xinjiang\s+Uighur\s+Autonomous\s+Region', loc_str, flags=re.IGNORECASE):
        loc_str = 'China,Xinjiang'
    elif re.search(r'Argentina', loc_str, flags=re.IGNORECASE) and re.search(r'Diamante.*Catamarca\s+Province', loc_str, flags=re.IGNORECASE):
        loc_str = 'Argentina,Catamarca'
    elif re.search(r'Chile', loc_str, flags=re.IGNORECASE) and re.search(r'South\s+America.*Lago\s+Ranco-Valdivia\s+Agricola\s+Quillin\s+Va\.\s+Region', loc_str, flags=re.IGNORECASE):
        loc_str = 'Chile' # not sure which location this should be
    elif re.search(r'Mexico', loc_str, flags=re.IGNORECASE) and re.search(r'Monarch\s+Butterfly\s+Biosphere\s+Reserve', loc_str, flags=re.IGNORECASE):
        loc_str = 'Mexico, Monarch Butterfly Biosphere Reserve'
    return loc_str

def preproc_loc_str(loc_str):
    """
    Pre-process location string
    Return processed result or None (for "no location")
    """
    import pandas
    if pandas.isnull(loc_str):
        return None
    # removes unicode or other "unexpected" chars
    loc_str = loc_str.encode('ascii','ignore').decode()
    # trailing/leading whitespaces
    loc_str = loc_str.strip()
    # strings known to encode "empty" entries
    if loc_is_missing(loc_str):
        return None
    # replace ":" by "," - helps in some cases to get coordinates
    loc_str = re.sub(':', ',', loc_str)
    # rm whitespaces after comma - helps in some cases to get coordinates
    loc_str = re.sub(',\s+', ',', loc_str)
    # exceptions
    loc_str = handle_loc_exceptions(loc_str)
    return loc_str

def add_dotzero(loc_str):
    """
    Add ".0" to a coordinate of it is only an int
    E.g. "23N" -> "23.0N"
    Required for correct query processing
    """
    if re.match(r'[0-9]+[N|S|W|E]', loc_str):
        loc_str = re.sub(r'([0-9]+)',r'\1.0',loc_str)
    return loc_str

def preproc_loc_coords(loc_str):
    """
    Pre-process location coordinate string
    Returns tuple of coordinates or None (for "no location")
    """
    import pandas
    if loc_str is None or pandas.isnull(loc_str):
        return None
    if re.fullmatch(r'[0-9]+(\.[0-9]+)?\s*[N|S][\s,]*[0-9]+(\.[0-9]+)?\s*[W|E]', loc_str):
        # matches = [add_dotzero(re.sub(' ','',m.group())) for m in re.finditer(r'[0-9]+(\.[0-9]+)?\s*[N|S|W|E]', loc_str)]
        matches = [re.sub(' ', '', m.group()) for m in re.finditer(r'[0-9]+(\.[0-9]+)?\s*[N|S|W|E]', loc_str)]
        assert len(matches) == 2, 'Problem to parse %s' % loc_str
        return matches[0], matches[1]
    elif re.fullmatch(r'[0-9]+(\.[0-9]+)?[\s,]*[0-9]+(\.[0-9]+)?', loc_str):
        matches = [re.sub(' ','',m.group()) for m in re.finditer(r'[0-9]+(\.[0-9]+)?', loc_str)]
        assert len(matches) == 2, 'Problem to parse %s' % loc_str
        return matches[0], matches[1]
    else:
        return None

def parse_location(loc_str, api_key, is_name=True):
    """
    Input: Location name or coordinates as string
    Returns {'lat': <float>, 'lng': <float>} if successful otherwise None
    Uses https://github.com/googlemaps/google-maps-services-python
    """
    import time, pandas, geocoder

    if pandas.isnull(loc_str):
        return {'lat': None, 'lng': None}

    if is_name:
        loc = geocoder.opencage(loc_str, key=api_key)
        time.sleep(1.5)
        # no hit
        if loc.latlng is None:
            logging.info('NO LOCATION HIT: %s' % loc_str)
            return {'lat': None, 'lng': None}
        else:
            return {'lat': loc.latlng[0], 'lng': loc.latlng[1]}
    else:
        assert len(loc_str) == 2
        lat, lng = loc_str
        lat_ = float(re.sub('N|S|W|E', '', lat))
        lng_ = float(re.sub('N|S|W|E', '', lng))
        if re.search('W', lng):
            lng_ *= -1
        if re.search('S', lat):
            lat_ *= -1
        return {'lat': lat_, 'lng': lng_}
