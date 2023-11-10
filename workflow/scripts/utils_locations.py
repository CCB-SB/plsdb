#!/usr/bin/python

import re
import logging
from shutil import ExecError
import pandas

##################################################
# LOCATIONS
##################################################
# Parsing location information from NCBI (BioSample)

def update_locs(locs, dict):
    import pandas as pd

    if locs is None:
        locs = pd.DataFrame([dict])
    else:
        locs = pd.concat([locs, pd.DataFrame([dict])])
    locs = locs.drop_duplicates(subset=["location"])
    locs.set_index('location', drop=False, inplace=True, verify_integrity=True)
    return locs

def save_locs(locs, ofile):
    assert locs is not None
    locs.to_csv(ofile, header=True, index=False, index_label=False, na_rep='')
    return


def loc_is_missing(loc_str):
    import re
    from utils import get_missing_info

    location_missing = get_missing_info()

    for na_str in location_missing:
        if re.fullmatch(na_str, loc_str, flags=re.IGNORECASE):
            return True
    return False

def handle_loc_exceptions(loc_str):
    import re
    if re.search(r'China,moutai town', loc_str, flags=re.IGNORECASE):
        return 'China,Maotaizhen'
    if re.search(r'Indonesia,East Java Province,Pacitan District,Donorojo Subdistrict,Klepu Village,Nglesen backwoods', loc_str):
        return 'Indonesia,East Java,Pacitan,Donorojo,Klepu,Nglesen'
    if re.search(r'China,near the Central Avenue of Tianjin Binhai', loc_str):
        return 'China,Binhai'
    if re.search(r"El Salvador,rural region", loc_str):
        return 'El Salvador'
    if re.search(r"Norway,Tautra underwater ridge at Trondheim Fjord", loc_str):
        return 'Norway,Trondheim Fjord'
    if re.search(r"China,the Shifeng River,Zhejiang Province", loc_str):
        return "China,Zhejiang,Taizhou,Shifeng Xi"
    if re.search(r"Norway,Troll Wall", loc_str):
        return "Norway,Trollveggen"
    if re.search(r"Chile,Microbial Pathogenesis and Vaccinology Laboratory,Huechuraba", loc_str):
        return "Chile,Huechuraba"
    if re.search(r"USA,Utqiagvik,Alaska CRREL Permafrost Tunnel", loc_str):
        return 'USA,AK99712,Fairbanks,Stesse Hwy,2126-2166'
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
    elif re.search('China,Zhejiang', loc_str, flags=re.IGNORECASE) and re.search('maternal and child health hospital', loc_str, flags=re.IGNORECASE):
        return 'China,Zhejiang,Maternal and Child Health Hospital'
    elif re.search('USA,Lost Angeles', loc_str, flags=re.IGNORECASE):
        return 'USA,Los Angeles'
    elif re.search('Germany,Gottingen', loc_str, flags=re.IGNORECASE):
        return 'Germany,Goettingen'
    return loc_str

def preproc_loc_str(loc_str):
    """
    Pre-process location string
    Return processed result or None (for "no location")
    """
    import re

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

# def add_dotzero(loc_str):
#     """
#     Add ".0" to a coordinate of it is only an int
#     E.g. "23N" -> "23.0N"
#     Required for correct query processing
#     """
#     if re.match(r'[0-9]+[N|S|W|E]', loc_str):
#         loc_str = re.sub(r'([0-9]+)',r'\1.0',loc_str)
#     return loc_str

def preproc_loc_coords(loc_str):
    """
    Pre-process location coordinate string
    Returns tuple of coordinates or None (for "no location")
    """
    import pandas
    import re

    if loc_str is None or loc_str == "nan":
        return None
    
    pattern = r'(\d+\.?\d+)\s*([NS])\s*;?(\d+\.?\d+)\s*([WE])'
    matches = re.search(pattern, loc_str)
    if matches:
        latitude = matches.group(1)
        latitude_direction = matches.group(2)
        longitude = matches.group(3)
        longitude_direction = matches.group(4)

        lat_str = latitude+latitude_direction
        lng_str = longitude+longitude_direction
        
        return [lat_str, lng_str]
    else:
        return None

def geopy_params_Nominatim(name = "location_mapping", language = "en", timeout=1):
    """
    https://geopy.readthedocs.io/en/stable/#specifying-parameters-once
    """
    from functools import partial
    from geopy.geocoders import Nominatim # Nominatim = OpenStreet Map
    
    geolocator = Nominatim(user_agent=name, timeout=timeout)
    geocode = partial(geolocator.geocode, language=language)
    
    
    return geocode

def geopy_params_Google(apikey, language = "en", timeout=1):
    """
    https://geopy.readthedocs.io/en/stable/#specifying-parameters-once
    """
    from functools import partial
    from geopy.geocoders import GoogleV3
    
    geolocator = GoogleV3(api_key=apikey, timeout=timeout)
    geocode = partial(geolocator.geocode, language=language)
    
    return geocode

def parse_location_geopy(geocode, loc_str, my_logger, is_name = True, country_code = True):
    """
    https://geopy.readthedocs.io/en/stable/#specifying-parameters-once
    """

    import pandas as pd
    import time

    if not loc_str or loc_str == "nan":
        return {'location': loc_str, 'lat': None, 'lng': None}
    
    if is_name:
        if country_code:
            # Enhance research by limiting results to the country
            list_str = loc_str.split(",")
            if len(list_str) > 1:
                country = list_str[0]
                country_code = get_country_code(country=country, 
                                                my_logger=my_logger)
            else:
                country_code = None

            if country_code:
                location = geocode(loc_str, country_codes = country_code)
            else:
                location = geocode(loc_str)
        else:
            location = geocode(loc_str)

        time.sleep(1.2)
        
        if location:
            my_logger.debug(f"{loc_str}: ({location.latitude} {location.longitude})")
            return {'location': loc_str,
                    'location_api': location.address,
                    'lat': location.latitude, 
                    'lng': location.longitude}
        else:
            return {'location': loc_str, 'lat': None, 'lng': None}
    else:
        if not len(loc_str) == 2:
            print(f"loc_str != 2 :{loc_str}")
            raise AssertionError
        lat, lng = loc_str
        lat_ = float(re.sub('N|S|W|E', '', lat))
        lng_ = float(re.sub('N|S|W|E', '', lng))
        if re.search('W', lng):
            lng_ *= -1
        if re.search('S', lat):
            lat_ *= -1
        return {'location': ';'.join(loc_str), 'lat': lat_, 'lng': lng_}
    
def interchange_positions(loc_str, my_logger):
    """
    Interchange second with third string. Sometimes helps to the search
    """ 
    loc_list = loc_str.split(',')

    if len(loc_list) >= 3:
        loc_elements = [loc_list[i] for i in [0,2,1]]

        if len(loc_list) > 3:
            loc_elements.extend(loc_list[3:len(loc_list)])
        
        # Join str
        loc_str_mod = ','.join(loc_elements)
        my_logger.info(f"Intercanging 2nd and 3rd string position: {loc_str} -> {loc_str_mod}")

        return loc_str_mod
    else:
        return None

def parse_location_google(geocode, loc_str):
    res = geocode(loc_str)

    if res:
        return {"location": loc_str,
                "location_api": res.address,
                "lat": res.latitude,
                "lng": res.longitude}

def parse_location_opencage(loc_str, api_key, is_name=True, logger=None):
    """
    Input: Location name or coordinates as string
    Returns {'lat': <float>, 'lng': <float>} if successful otherwise None
    Uses https://github.com/googlemaps/google-maps-services-python
    """
    import time, pandas, geocoder

    if pandas.isnull(loc_str) or loc_str=="nan":
        return {'location': loc_str, 'lat': None, 'lng': None}

    if is_name:
        loc = geocoder.opencage(loc_str, key=api_key)
        time.sleep(1.5)
        # no hit
        if loc.latlng is None:
            # logger.info(f'NO LOCATION HIT: {loc_str}')
            return {'location': loc_str, 'lat': None, 'lng': None}
        else:
            return {'location': loc_str,
            'lat': loc.latlng[0], 'lng': loc.latlng[1]}
    else:
        if not len(loc_str) == 2:
            print(f"loc_str != 2 :{loc_str}")
            raise AssertionError
        lat, lng = loc_str
        lat_ = float(re.sub('N|S|W|E', '', lat))
        lng_ = float(re.sub('N|S|W|E', '', lng))
        if re.search('W', lng):
            lng_ *= -1
        if re.search('S', lat):
            lat_ *= -1
        return {'location': loc_str, 'lat': lat_, 'lng': lng_}

def fuzzy_search_locations(loc_str, loc_df,  my_logger, threshold=80):
    """
    Make a fuzzy search against our location database
    """
    from rapidfuzz import process

    # If existing mapping
    # If match in existing mapping with >= 95% similarity
    match, score, choice = process.extractOne(loc_str, list(loc_df['location']))
    
    if match and score >= threshold:
        res = {
            "location": match,
            'lat': loc_df.loc[match,'lat'],
            'lng': loc_df.loc[match,'lng']}
        my_logger.debug(f"Fuzzy match: {loc_str} -> {match}")
    else:
        res = {'location': loc_str, 'lat': None, 'lng': None}
    
    return res

def get_country_code(country, my_logger):
    import pycountry

    try:
        search = pycountry.countries.search_fuzzy(country)
        code = search[0].alpha_2.lower()
        return code
    except Exception as e:
        my_logger.warning(f"No country code for: {country}")
        return None

def compare_gecoder_results(nominatim, google):
    from geopy.distance import geodesic as GD
    from rapidfuzz import fuzz

    # Compare strings:
    s1 = nominatim["location_api"]
    s2 = google["location_api"]
    ratio = fuzz.token_set_ratio(s1, s2)

    # Compare locations:
    l1 = (nominatim["lat"], nominatim['lng'])
    l2 =(google["lat"] , google['lng'])
    dist = GD(l1,l2).km

    return ratio, dist

def compute_middle_point(locs):
    """
    Compute middle point between locations. Useful when more than one
        city is specified.
    
    ::params locs: list of str
    ::return (lat, lng) of middle point
    """
    from utils_my_logger import My_logger
    log = My_logger(log_filename="tmp.log").get_logger()
    
    geocode_nom = geopy_params_Nominatim(name="middle_point", language="en", timeout=2)
    
    lat = []
    lng = []
    for loc in locs:
        res = parse_location_geopy(geocode_nom, loc, is_name=True, my_logger=log, 
                                   country_code=False)
        lat.append(res['lat'])
        lng.append(res['lng'])
    
    mid_lat = sum(lat)/len(lat)
    mid_lng = sum(lng)/len(lng)

    return (mid_lat, mid_lng)

def reverse_string(s):
    s_list = s.split(',')

    if len(s_list) > 1:
        s_list.reverse()
        s = ','.join(s_list)

    return s
    

def reverse_output_location(df, version_name):
    """
    Reverse order of location address.
    Useful to only copy paste locations from browser and 
        perform reversion name afterwards.

    ::params df (input a copy of the df to avoid problems)       
    """
    import numpy as np

    df['reverse'] = np.where((df['version'] == version_name) & (df['notes'].str.contains('coordinates')),
              False,
              np.where((df['version'] == version_name),
                       True,
                       False))
    
    df.loc[df['reverse']==True,'output_location'] = df.loc[df['reverse']==True, 'output_location'].apply(reverse_string)

    return df['output_location']

