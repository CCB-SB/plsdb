
from ast import parse
from re import A
import utils_locations
from utils_my_logger import My_logger
from utils import unify_missing_info, create_mapping
import pandas as pd
import numpy as np

##################################################
# ARGS
##################################################
def get_arg_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--plasmids', '-p', help="File with all plasmids.", required=True)
    parser.add_argument('--ofile', '-o', help="Output file for locations.", required=True)
    parser.add_argument('--output-mapping', help="Output file for new locations.", required=True)
    parser.add_argument('--output-nohits', help="Output file with locations without hits", required=True)
    parser.add_argument("--google-api", required = True, help="Google API KEY")
    parser.add_argument("--version", required = True, help="new version str")
    parser.add_argument("--log", required = True)
    parser.add_argument("--locs-correction", required = True,  help = "File to correct input location names")
    parser.add_argument('--locs', '-l', required=True,
                        help="Path to file for (already retrieved) locations.")
    
    return parser

##################################################
# FUNCTIONS
##################################################

##################################################
# MAIN
##################################################

if __name__ == "__main__":
    ARGS = get_arg_parser().parse_args()
    # Logging
    logger = My_logger(log_filename = ARGS.log, logger_name = "infer_location")
    my_logger = logger.get_logger()

    # API key
    api_key = ARGS.google_api

    # load known locations
    loc_df = pd.read_csv(ARGS.locs,na_values=[''])
    loc_df.set_index('location', drop=False, inplace=True, verify_integrity=True)
    loc_df.dropna(subset=["lat", "lng"], inplace=True)
    my_logger.info(f'Loaded locations from {ARGS.locs}\n{loc_df.head()}')

    # Load correction for locations
    loc_corrections = pd.read_csv(ARGS.locs_correction, dtype=str, na_values=[''])
    correction_dict = create_mapping(loc_corrections, 
                   input_column="input_location", 
                   output_column="output_location")
    ## check if the corrected location is directly a coordinate
    correction_coor= create_mapping(loc_corrections, 
                   input_column="input_location", 
                   output_column="notes")
    
    # Configure geopy
    geocode_nom = utils_locations.geopy_params_Nominatim(name="location_mapping",
                                           language="en", timeout=2)
    geocode_go = utils_locations.geopy_params_Google(apikey=ARGS.google_api,
                                                  language='en', timeout=2)

    # sample table with locations
    pls = pd.read_csv(ARGS.plasmids, dtype=str)
    pls.set_index('NUCCORE_UID', drop=False, inplace=True, verify_integrity=True)
    pls = unify_missing_info(pls, value_missing="nan", my_logger=my_logger)

    # parse
    coords = []
    no_hits = []

    # For location name and coordinates
    for ID, loc, coor in zip(list(pls["NUCCORE_UID"]), 
                             list(pls["BIOSAMPLE_Location"]), 
                             list(pls["BIOSAMPLE_Coordinates"])):

        # no location data
        if (loc is None or loc == "nan") and (coor is None or coor == "nan"):
            continue
        
        # pre-processing
        # my_logger.debug("Preprocessing loc and coor...")
        # my_logger.debug(f"loc: '{loc}'")
        loc = utils_locations.preproc_loc_str(loc)
        # my_logger.debug(f"loc preprocessed: '{loc}'")
        # my_logger.debug(f"coor: {coor}")
        coor = utils_locations.preproc_loc_coords(utils_locations.preproc_loc_str(coor))
        # my_logger.debug(f"coor preprocessed: {coor}")
        
        if (loc is None or loc == "nan") and (coor is None or coor == "nan"):
            continue
        
        if loc in correction_dict.keys():
            if 'coordinates' in str(correction_coor[loc]):
                coor = correction_dict[loc].split(';')
                my_logger.debug(f"Adding coordinates from corrections: loc = '{coor}'")
            else:
                loc = correction_dict[loc]
                my_logger.debug(f"Adding location from corrections: loc = '{loc}'")
        
        try:
            # at least name or coordinates given
            # If coordinate
            if coor is not None and coor != "nan":
                coor_str = f'{coor[0]};{coor[1]}'
                if coor_str not in loc_df['location']:
                    my_logger.info(f'Retrieving coordinates for location: {coor_str}')

                    # Search in OpenStreet Map
                    parsed = utils_locations.parse_location_geopy(
                        geocode_nom, loc_str=coor, is_name=False, my_logger=my_logger)
                    
                    # Update local db if result
                    if parsed['lat'] and parsed['lng']:
                        loc_df = utils_locations.update_locs(loc_df, 
                                        {'location': parsed['location'], 'type': 'coordinates', 
                                            'lat': parsed['lat'], 'lng': parsed['lng'],
                                            'version': ARGS.version})
                        utils_locations.save_locs(loc_df, ARGS.output_mapping)
                        # Result for plasmid table
                        coords.append({'ID': ID, 
                                        'loc_lat': parsed['lat'], 
                                        'loc_lng': parsed['lng']})
                    else:
                        my_logger.info(f'NO LOCATION HIT: {coor_str}')
                        no_hits.append({
                            'input_location': coor_str, 
                            'output_location': "",
                            "version": ARGS.version})
                    
                else:
                    coords.append({'ID': ID, 'loc_lat': loc_df.loc[coor_str,'lat'], 
                                'loc_lng': loc_df.loc[coor_str,'lng']})

            # Or location name
            elif loc is not None and loc != "nan":
                # my_logger.debug(f"Empty coor: {coor}")
                if loc not in loc_df['location']:
                    fuzzy = False
                    my_logger.info(f'Retrieving coordinates for location: \"{loc}\"')


                    # Search in OpenStreet Map
                    parsed = utils_locations.parse_location_geopy(geocode_nom,
                        loc_str = loc, is_name=True, my_logger=my_logger)
                    
                    # Fuzzy search against local db: prevent typos
                    if not parsed['lat'] and not parsed['lng']:
                        parsed = utils_locations.fuzzy_search_locations(
                            loc_str=loc, loc_df = loc_df, threshold = 95,
                            my_logger = my_logger)
                        fuzzy = True
                    
                    # Invert locations
                    if not parsed['lat'] and not parsed['lng']:
                        loc_mod = utils_locations.interchange_positions(loc_str=loc, 
                                                                        my_logger=my_logger)
                        if loc_mod:
                            # Fuzzy search
                            parsed = utils_locations.fuzzy_search_locations(
                                loc_str=loc_mod, loc_df = loc_df, threshold = 95,
                                my_logger = my_logger)
                            fuzzy = False
                            
                            if not parsed['lat'] and not parsed['lng']:
                                # Search OpenStreet
                                parsed = utils_locations.parse_location_geopy(geocode_nom,
                                    loc_str = loc_mod, is_name=True, my_logger=my_logger)
                                fuzzy = True

                    # if result
                    if parsed['lat'] and parsed['lng']:
                        if not fuzzy:
                            # my_logger.debug(f"fuzzy: {fuzzy}")
                            # Search in Google to compare
                            parsed_go = utils_locations.parse_location_google(geocode_go,
                                loc_str = loc)
                            
                            if parsed_go:
                                my_logger.debug("Simarility Nominatim vs Google:")
                                sim, dist = utils_locations.compare_gecoder_results(
                                    nominatim=parsed, google=parsed_go)
                                my_logger.debug(f"String locations = {sim}")
                                my_logger.debug(f"Distance = {dist} km")
                            else: 
                                sim, dist = None, None
                            
                            # Update local db
                            loc_df = utils_locations.update_locs(loc_df, 
                                                        {'location': parsed['location'], 
                                                        'location_api': parsed['location_api'],
                                                        'dist_km': dist, 'similarity_locs': sim, 
                                                        'type': 'name', 
                                                        'lat': parsed['lat'], 'lng': parsed['lng'],
                                                        'version': ARGS.version})
                            utils_locations.save_locs(loc_df, ARGS.output_mapping)
                        
                        # Result for plasmid table
                        coords.append({'ID': ID, 'loc_lat': parsed['lat'], 
                                    'loc_lng': parsed['lng'],
                                    'loc_parsed': parsed['location']})
                    else:
                        my_logger.info(f'NO LOCATION HIT: {loc}')
                        no_hits.append({
                            'input_location': loc, 
                            'output_location': "",
                            "version": ARGS.version})
                    
                else:
                    coords.append({'ID': ID, 'loc_lat': loc_df.loc[loc,'lat'], 
                                'loc_lng': loc_df.loc[loc,'lng'],
                                'loc_parsed': loc})
    
        except Exception as e:
            # Save no hits information
            nohits_df = pd.DataFrame(no_hits)
            nohits_df.drop_duplicates(inplace=True)
            nohits_df = pd.concat([loc_corrections, nohits_df])
            nohits_df.to_csv(ARGS.output_nohits, index=False, quoting=2)
            raise e
    
    # Save plasmid table
    coords = pd.DataFrame(coords)
    coords.set_index('ID', drop=True, verify_integrity=True, inplace=True)
    pls = pd.merge(left=pls, right=coords, how='left',
                   left_index=True, right_index=True, sort=False)

    pls.to_csv(ARGS.ofile, header=True, index=False)

    # Save no hits information
    nohits_df = pd.DataFrame(no_hits)
    nohits_df.drop_duplicates(inplace=True)
    nohits_df = pd.concat([loc_corrections, nohits_df])
    nohits_df.to_csv(ARGS.output_nohits, index=False, quoting=2)
