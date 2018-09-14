#!/usr/bin/python

import re
import logging
import pandas

##################################################
# MISC
##################################################
def setup_logger(log_level=logging.DEBUG, log_file='pipeline.log'):
    """
    Set up the logger
    """
    # file handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(log_level)
    # console handler
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    # formatter
    formatter = logging.Formatter(
        fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%y:%m:%d %H:%M:%S'
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # logger
    logger = logging.getLogger('pipeline_logger')
    logger.setLevel(log_level)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

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
    logger = setup_logger(logging.INFO)

    logger.info('CMD: {}'.format(cmd))

    ids = pandas.read_csv(df_file, sep='\t', header=0, dtype=str)[df_col]
    logger.info('There are {} IDs'.format(len(ids)))
    ids.fillna(value='', inplace=True)
    ids = ids.loc[ids != ""]
    ids = set(ids.values)
    if split_str is not None and split_str != "":
        ids_ = set()
        for id_ in ids:
            ids_ = ids_.union(id_.split(split_str))
        ids = ids_
    logger.info('There are {} unique IDs'.format(len(ids)))

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
                    logger.warn('Could not find: %s' % ', '.join([str(s) for s in not_found]))
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

    logger = setup_logger(logging.INFO)

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
        header=' '.join(blastn_header),
        pident=blastn_pident
    )
    logger.info('START {}\n{}'.format(acc, cmd))
    cmd, cmd_s, cmd_o = run_cmd(cmd)
    assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd, cmd_s, cmd_o)
    logger.info('END {}'.format(acc))
    return

def proc_mlst_scheme_name(scheme_name):
    """
    Process scheme name, e.g. to be used as directory name
    """
    import re
    return re.sub('/', '_', re.sub('\s+', '__', scheme_name))

def reproc_mlst_scheme_name(scheme_name):
    """
    Reverts the changes made by proc_mlst_scheme_name()
    """
    import re
    return re.sub('_', '/', re.sub('__', ' ', scheme_name))

def download_pmlst_scheme_alleles(scheme_name, scheme_url, scheme_dir):
    """
    Download pMLST scheme allele sequences
    :param scheme_name: name of the scheme
    :param: scheme_url: URL to be used
    :param scheme_dir: where to save the output files
    """
    import os
    import requests

    logger = setup_logger(logging.INFO)

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
            logger.warning('Scheme {}, locus {}: no alleles'.format(scheme_name, locus))
    return

def download_pmlst_scheme_profiles(scheme_name, scheme_url, scheme_dir, scheme_profiles):
    """
    Download pMLST scheme profiles for given scheme
    Will also check the formatting s.t. the files can be used by the "mlst" tool
    :param scheme name: name of the scheme
    :param scheme_url: URL to be used
    :param scheme_dir: where the allele FASTA files (*.tfa) are stored
    :param scheme_profiles: output file name for the profiles
    """
    import os
    import re
    import requests
    from glob import glob

    logger = setup_logger(logging.INFO)

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

##################################################
# LOCATIONS
##################################################
# Parsing location information from NCBI (BioSample)

def load_locs(loc_file):
    import pandas, os

    logger = setup_logger(logging.INFO)

    locs = None
    if os.path.exists(loc_file):
        locs = pandas.read_csv(loc_file, sep='\t', header=0)
        locs.set_index('location', drop=False, inplace=True, verify_integrity=True)
        logger.info('Loaded locations from {}\n{}'.format(loc_file, locs.head()))
    else:
        logger.info('No known locations to load.')
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
    if re.search(r'BIFSCo Region \d?', loc_str, flags=re.IGNORECASE) and re.search(r'US|USA|United States', loc_str, flags=re.IGNORECASE):
        return 'USA'
    elif re.search(r'^Australia,Sealake$', loc_str, flags=re.IGNORECASE):
        return re.sub('sealake', 'Sea Lake', loc_str)
    elif re.search(r'^South Korea,Yeo-su sediment$', loc_str, flags=re.IGNORECASE):
        return 'South Korea,Yeosu'
    elif re.search(r'^South Korea,Shnan-gun$', loc_str, flags=re.IGNORECASE):
        return 'South Korea'
    elif re.search(r'South Korea,Geojedo', loc_str, flags=re.IGNORECASE):
        return 'South Korea,Geoje'
    elif re.search(r'South Korea,the surface of the seashore around a seaweed farm at Geoje Island in the South Sea', loc_str, flags=re.IGNORECASE):
        return 'South Korea,Geoje'
    elif re.search(r'^China,Harbin Veterinary Research Institute$', loc_str, flags=re.IGNORECASE):
        return 'China,Harbin'
    elif re.search(r'Germany', loc_str, flags=re.IGNORECASE) and re.search(r'Black Forest', loc_str, flags=re.IGNORECASE):
        return 'Germany,Black Forest'
    elif re.search(r'^Sweden,Kosterfjord$', loc_str, flags=re.IGNORECASE):
        return 'Sweden,Kosterfjorden'
    elif re.search(r'^China,Bo Hai,Panjin$', loc_str, flags=re.IGNORECASE):
        return 'China,Panjin'
    elif re.search(r'^China,Sunitezuoqi$', loc_str, flags=re.IGNORECASE):
        return 'China' # location unknown even in Google Maps
    elif re.search(r'^Brazil,State of Rondonia,Western Amazon$', loc_str, flags=re.IGNORECASE):
        return 'Brazil,Rondonia'
    elif re.search(r'^South Korea,Tae-an sediment$', loc_str, flags=re.IGNORECASE):
        return 'South Korea,Taean'
    elif re.search(r'^South Africa,Kruger National Park,Pafuri$', loc_str, flags=re.IGNORECASE):
        return 'South Africa,Kruger National Park'
    elif re.search(r'^Japan,Hyogo,Himeji,Univrersity of Hyogo,Harima Campus for Science$', loc_str, flags=re.IGNORECASE):
        return 'Japan,Himeji'
    elif re.search(r'^USA,Texas,Baytown,Burnet Shores$', loc_str, flags=re.IGNORECASE):
        return 'USA,Baytown'
    elif re.search(r'^Thailand,Nakhornrachisma$', loc_str, flags=re.IGNORECASE):
        return 'Thailand,Nakhon Ratchasima'
    elif re.search(r'^China,Eastern Hubei Province$', loc_str, flags=re.IGNORECASE):
        return 'China,Hubei'
    elif re.search(r'^China,Chenmai qiaotou$', loc_str, flags=re.IGNORECASE):
        return 'China,Qiaotouzhen'
    elif re.search(r'^China,Hongyuan Prairie,Aba Autonomous Prefecture Homo sapiens$', loc_str, flags=re.IGNORECASE):
        return 'China,Hongyuan'
    elif re.search(r'^Brazil,Tupasi$', loc_str, flags=re.IGNORECASE):
        return 'Brazil'
    elif re.search(r'China', loc_str, flags=re.IGNORECASE) and re.search(r'Xinjiang Province', loc_str, flags=re.IGNORECASE):
        return 'China,Xinjiang'
    elif re.search(r'^USA,Utah,Bear River Refuge$', loc_str, flags=re.IGNORECASE):
        return 'USA,Bear River Migratory Bird Refuge' # the only matching refuge location
    elif re.search(r'^USA,Utah,Ogden Bay Refuge$', loc_str, flags=re.IGNORECASE):
        return 'USA,Ogden'
    elif re.search(r'^Vietnam,Do Xongpha$', loc_str, flags=re.IGNORECASE):
        return 'Vietnam' # could not find the associated location
    elif re.search(r'^India,Haiderabad/Hindustan$', loc_str, flags=re.IGNORECASE):
        return 'India,Haiderabad'
    elif re.search(r'^Democratic Republic of the Congo,Iturie province$', loc_str, flags=re.IGNORECASE):
        return 'Democratic Republic of the Congo,Ituri'
    elif re.search(r'^Germany,Bruschal$', loc_str, flags=re.IGNORECASE):
        return 'Germany,Bruchsal' # typo?
    elif re.search(r'^Canada,University of British Columbia,pilot reactor for wastewater decontamination$', loc_str, flags=re.IGNORECASE):
        return 'Canada,Vancouver'
    elif re.search(r'^Brazil,Ribeirao Preto,Sao Paulo State$', loc_str, flags=re.IGNORECASE):
        return 'Brazil,Ribeirao Preto'
    elif re.search(r'^Japan,Niigata,Nagakura$', loc_str, flags=re.IGNORECASE):
        return 'Japan,Nagakura'
    elif re.search(r'^China,Xinjiang Uighur Autonomous Region$', loc_str, flags=re.IGNORECASE):
        return 'China,Xinjiang'
    elif re.search(r'^Argentina,Diamante,Catamarca Province$', loc_str, flags=re.IGNORECASE):
        return 'Argentina,Catamarca'
    elif re.search(r'^Chile,South America,Lago Ranco-Valdivia Agricola Quillin Va\. Region$', loc_str, flags=re.IGNORECASE):
        return 'Chile' # not sure which location this should be
    elif re.search(r'^Mexico,Monarch Butterfly Biosphere Reserve,Sierra Chivati-Huacal$', loc_str, flags=re.IGNORECASE):
        return 'Mexico,Monarch Butterfly Biosphere Reserve'
    elif re.search(r'^Antarctica,1670 km from geographic south pole$', loc_str, flags=re.IGNORECASE):
        return 'Antarctica'
    elif re.search(r'^Belarus,Minsk area,Myadzel reg\.$', loc_str, flags=re.IGNORECASE):
        return 'Belarus,Myadzel'
    elif re.search(r'^Chile,X Region$', loc_str, flags=re.IGNORECASE):
        return 'Chile,Region de los Lagos'
    elif re.search(r'^China,coastal sediment close to a coal-fired power station,Qingdao$', loc_str, flags=re.IGNORECASE):
        return 'China,Qingdao'
    elif re.search(r'^China,General Microbiological Culture Collection Center$', loc_str, flags=re.IGNORECASE):
        return 'China,Beijing'
    elif re.search(r'^China,Jilin Pesticide Plant$', loc_str, flags=re.IGNORECASE):
        return 'China,Jilin'
    elif re.search(r'^China,Jilin,Jilin Oil Field Branch$', loc_str, flags=re.IGNORECASE):
        return 'China,Jilin'
    elif re.search(r'^China,western desert$', loc_str, flags=re.IGNORECASE):
        return 'China'
    elif re.search(r'^Germany,cities of Ulm and Neu-Ulm,south-west Germany$', loc_str, flags=re.IGNORECASE):
        return 'Germany,Ulm,Neu-Ulm'
    elif re.search(r'^Guangxi,China$', loc_str, flags=re.IGNORECASE):
        return 'China,Guangxi'
    elif re.search(r'^Indonesia,hot spring$', loc_str, flags=re.IGNORECASE):
        return 'Indonesia'
    elif re.search(r'^Jamaica,Hector\'s Bay$', loc_str, flags=re.IGNORECASE):
        return 'Jamaica'
    elif re.search(r'^Japan,Kagoshima,Tokunoshima Island,Rikuhama beach$', loc_str, flags=re.IGNORECASE):
        return 'Japan,Kagoshima,Tokunoshima Island'
    elif re.search(r'^korea,chung-nam$', loc_str, flags=re.IGNORECASE):
        return 'Korea,Chungnam'
    elif re.search(r'^Mexico,Huautla,San Miguel Acuexcomac,Pueblo$', loc_str, flags=re.IGNORECASE):
        return 'Mexico,Pueblo,San Miguel Acuexcomac'
    elif re.search(r'^Netherlands,Bennekom area$', loc_str, flags=re.IGNORECASE):
        return 'Netherlands,Bennekom'
    elif re.search(r'^Norway,Skammestein at route 51 northwest of Fagernes,intersection to E16$', loc_str, flags=re.IGNORECASE):
        return 'Norway,Fagernes'
    elif re.search(r'^Oman,Jabal Al-Akhdar$', loc_str, flags=re.IGNORECASE):
        return 'Oman,Al Jabal Al Akhdar'
    elif re.search(r'^Russia,Baikal Lake region,Siberia$', loc_str, flags=re.IGNORECASE):
        return 'Russia,Baikal Lake'
    elif re.search(r'^Russia,Kolyma,Lowlands,Siberia$', loc_str, flags=re.IGNORECASE):
        return 'Russia,Kolyma,Lowlands'
    elif re.search(r'^Russia,Nadym region of northwest Siberia,Yamalo-Nenets AO$', loc_str, flags=re.IGNORECASE):
        return 'Russia,Nadym'
    elif re.search(r'^Russia,Pacific Ocean,Rudnaya Bay$', loc_str, flags=re.IGNORECASE):
        return 'Russia,Bukhta Rudnaya'
    elif re.search(r'^South Africa$', loc_str, flags=re.IGNORECASE):
        return 'South Africa,Africa'
    elif re.search(r'^South Korea,Cheju Island$', loc_str, flags=re.IGNORECASE):
        return 'South Korea,Jeju-do'
    elif re.search(r'^South Korea,chung-nam$', loc_str, flags=re.IGNORECASE):
        return 'South Korea,chungnam'
    elif re.search(r'^Spain,Malaga,Estacion Experimental "La Mayora",Algarrobo Costa$', loc_str, flags=re.IGNORECASE):
        return 'Spain,Algarrobo Costa'
    elif re.search(r'^Taiwan,Kueishan Island$', loc_str, flags=re.IGNORECASE):
        return 'Taiwan,Kueishan'
    elif re.search(r'^Ukraine,Crimean Lake Mainaki silt$', loc_str, flags=re.IGNORECASE):
        return 'Ukraine,Crimea'
    elif re.search(r'^USA,New York ex Aruba$', loc_str, flags=re.IGNORECASE):
        return 'USA,New York'
    elif re.search(r'^USSR$', loc_str, flags=re.IGNORECASE):
        return None
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
    # replace ":" by ","
    loc_str = re.sub(':', ',', loc_str)
    # rm trailing comma
    loc_str = re.sub(',$', '', loc_str)
    # rm whitespaces after comma
    loc_str = re.sub(',\s+', ',', loc_str)
    # replace multiple whitespaces by one
    loc_str = re.sub('\s\s+', ' ', loc_str)

    # strings known to encode "empty" entries
    if loc_is_missing(loc_str):
        return None

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

    logger = setup_logger(logging.INFO)

    if pandas.isnull(loc_str):
        return {'lat': None, 'lng': None}

    if is_name:
        loc = geocoder.opencage(loc_str, key=api_key)
        time.sleep(1.5)
        # no hit
        if loc.latlng is None:
            logger.info('NO LOCATION HIT: %s' % loc_str)
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
