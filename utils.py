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
    import os
    if p:
        if not os.path.exists(dpath):
            os.makedirs(dpath)
    else:
        os.makedirs(dpath)

def proc_mlst_scheme_name(scheme_name):
    import re
    return re.sub('/', '_', re.sub('\s+', '__', scheme_name))

def reproc_mlst_scheme_name(scheme_name):
    import re
    return re.sub('_', '/', re.sub('__', ' ', scheme_name))

##################################################
# Edirect/eutils
##################################################
def run_epost_split(df_file, ofile, header, cmd, df_col, split_size, split_str, **kwargs):
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

##################################################
# rMLST
##################################################
def run_rmlst(sfasta):
    """
    Run rMLST analysis
    Reference: https://pubmlst.org/rmlst/api.shtml
    Ribosomal Multilocus Sequence Typing (rMLST) is an approach that indexes
    variation of the 53 genes encoding the bacterial ribosome protein subunits
    (rps genes) as a means of integrating microbial taxonomy and typing.
    The rps gene variation catalogued in this database permits rapid and
    computationally non-intensive identification of the phylogenetic position
    of any bacterial sequence at the domain, phylum, class, order, family, genus,
    species and strain levels.
    Reference: https://www.ncbi.nlm.nih.gov/pubmed/22282518
    The rps loci are [...]
    (i)   present in all bacteria;
    (ii)  distributed around the chromosome;
    (iii) encode proteins which are under stabilizing selection for functional conservation

    :param sfasta: FASTA file <path>/<seqID>.<ext>
    :return dictionary with accession ID, and rMLST support/rank/taxon/lineage (best hit w.r.t. support)
    """
    import os, sys, requests, base64
    from random import randint
    from time import sleep

    # check that FASTA file exists
    assert os.path.exists(sfasta), 'File {} does not exist'.format(sfasta)

    # derive sequence ID/accession
    sid = os.path.splitext(os.path.basename(sfasta))[0]

    # read in sequence data
    with open(sfasta, 'r') as ifile:
        payload = '{"base64":true,"details":true,"sequence":"' + base64.b64encode(ifile.read().encode()).decode() + '"}'
    # submit for analysis
    response = requests.post(
        'http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence',
        data=payload
    )
    sleep(1 + 1/randint(2, 10)) # to avoid too many queries per second
    # get best match
    best_match = {
        'ACC_FASTA': sid,
        'Rank_RMLST': None,
        'Taxon_RMLST': None,
        'Support_RMLST': None,
        'Taxonomy_RMLST': None,
    }
    assert response.status_code == requests.codes.ok, 'Response status not OKAY for {}: {}'.format(sid, response.status_code)
    try:
        for match in response.json()['taxon_prediction']:
            match = {
                'ACC_FASTA': sid,
                'Rank_RMLST': match['rank'],
                'Taxon_RMLST': match['taxon'],
                'Support_RMLST': match['support'],
                'Taxonomy_RMLST': match['taxonomy'],
            }
            if best_match['Support_RMLST'] is None or match['Support_RMLST'] > best_match['Support_RMLST']:
                best_match = match
        return best_match
    except KeyError:
        return best_match

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

# Exceptions: need to change string to get result from GMaps
location_exceptions = {
    # 'Jigalong, Australia': 'Australia, Jigalong', # ordering matters
    # 'China: Jilin Pesticide Plant': 'China: Jilin', # non-loc. info.
    # 'Korea: Gwangyang Province': 'Korea: Gwangyang', #
    # 'China:Pearl Spring': 'China', # unknown location
    # 'China: the China South Sea': 'South China Sea', # otherwise maps to US
    # 'China: Changji city, Xinjiang Uygur Autonomous Region': 'China: Changji', # maps to loc. in China near border to S: Korea
    # 'China: Southern China': 'China', # otherwise maps to US
}

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
        if re.fullmatch(na_str, loc_str, re.IGNORECASE):
            return True
    return False

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
    # exceptions
    if loc_str in location_exceptions:
        loc_str = location_exceptions[loc_str]
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

# def parse_location(loc_str, api_key, is_name=True):
#     """
#     Input: Location name or coordinates as string
#     Returns {'lat': <float>, 'lng': <float>} if successful otherwise None
#     Uses https://github.com/googlemaps/google-maps-services-python
#     """
#     import time
#     import googlemaps
#
#     if pandas.isnull(loc_str):
#         return None
#
#     input = loc_str
#     gmaps = googlemaps.Client(key=api_key)
#     loc = None
#     # query
#     loc = None
#     if is_name or re.search(r'N|W|E|S', loc_str[0]) or re.search(r'N|W|E|S', loc_str[1]):
#         loc = gmaps.geocode(loc_str)
#         time.sleep(1)
#     else:
#         loc = gmaps.reverse_geocode(loc_str)
#         time.sleep(1)
#     # no hit
#     if len(loc) == 0:
#         logging.info('Location: no hit: %s' % input)
#         return None
#     # got location
#     assert 'geometry' in loc[0], 'Location: hit w/o geometry: %s' % input
#     assert 'location' in loc[0]['geometry'], 'Location: hit w/o coordinates: %s' % input
#     return loc[0]['geometry']['location']
