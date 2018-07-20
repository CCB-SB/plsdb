#!/usr/bin/python

import os
import re
import logging

##################################################
# MISC
##################################################
def mkdir_p(p):
    import os
    if not os.path.exists(p):
        os.makedirs(p)

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

def setup_logger(level=logging.INFO):
    """Logger"""
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%y:%m:%d %H:%M:%S'
    )

def rm_file(fpath):
    if os.path.exists(fpath):
        os.remove(fpath)

def is_gzipped(fname):
    """
    Naive way to check whether a file is gzipped
    """
    return re.search('\.gz$', fname) or re.search('\.gzip', fname)

def chunk_list(ilist, num=1, length=None):
    """
    Split list into chunks
    """
    import math
    ilist = list(ilist)
    assert num is not None or length is not None
    if num is not None:
        length = math.ceil(len(ilist) / num)
    for i in range(0, len(ilist), length):
        yield ilist[i:i + length]

##################################################
# LOCATIONS
##################################################
# Parsing location information from NCBI (BioSample) using GoogleMaps

# Exceptions: need to change string to get result from GMaps
location_exceptions = {
    'Jigalong, Australia': 'Australia, Jigalong',
    'China: Jilin Pesticide Plant': 'China: Jilin',
    'Korea: Gwangyang Province': 'Korea: Gwangyang',
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

def preproc_loc_str(loc_str):
    """
    Pre-process location string
    Return processed result or None (for "no location")
    """
    import pandas
    if loc_str is None or pandas.isnull(loc_str):
        return None
    # removes unicode or other "unexpected" chars
    loc_str = loc_str.encode('ascii','ignore').decode()
    # trailing/leading whitespaces
    loc_str = loc_str.strip()
    # strings known to encode "empty" entries
    for na_str in location_missing:
        if re.fullmatch(na_str, loc_str, re.IGNORECASE):
            return None
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
        matches = [add_dotzero(re.sub(' ','',m.group())) for m in re.finditer(r'[0-9]+(\.[0-9]+)?\s*[N|S|W|E]', loc_str)]
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
    import time
    import googlemaps
    input = loc_str
    gmaps = googlemaps.Client(key=api_key)
    loc = None
    # exceptions
    if loc_str in location_exceptions:
        loc_str = location_exceptions[loc_str]
    # pre-processing
    loc_str = preproc_loc_str(loc_str)
    if not is_name:
        loc_str = preproc_loc_coords(loc_str)
    # empty
    if loc_str is None:
        logging.info('Location: no hit: %s' % input)
        return None
    # query
    loc = None
    if is_name or re.search(r'N|W|E|S', loc_str[0]) or re.search(r'N|W|E|S', loc_str[1]):
        loc = gmaps.geocode(loc_str)
        time.sleep(1)
    else:
        loc = gmaps.reverse_geocode(loc_str)
        time.sleep(1)
    # no hit
    if len(loc) == 0:
        logging.info('Location: no hit: %s' % input)
        return None
    # got location
    assert 'geometry' in loc[0], 'Location: hit w/o geometry: %s' % input
    assert 'location' in loc[0]['geometry'], 'Location: hit w/o coordinates: %s' % input
    return loc[0]['geometry']['location']
