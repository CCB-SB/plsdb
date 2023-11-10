#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 
# Authors:   Alejandra Gonzalez (gonmola@hotmail.es)
# ---

##################################################
# IMPORTS
##################################################
from utils_my_logger import My_logger
from utils import create_mapping, unify_missing_info
import pandas as pd
import os
##################################################
# ARGS
##################################################


def get_parser():
    import argparse
    import json
    # create parser object
    parser = argparse.ArgumentParser(description="Create dictionary like tsv files to facilitate the host inferring process",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # specify our desired options 
    # by default ArgumentParser will add an help option 
    parser.add_argument("--input-pls", help="CSV metadata file of plasmids")
    parser.add_argument("--host-mapping", help="csv file for host mapping from BIOSAMPLE_Host")
    parser.add_argument("--outfile-host-mapping", help="csv outfile name for host mapping")
    parser.add_argument("--outfile-host-mapping-checked", help="csv outfile name for host mapping")
    parser.add_argument("--version", 
                        help="Name for the newversion of the mapping files. Please follow the format Y/m/d (e.g 2023_11_01)")
    parser.add_argument("--api-key", help="NCBI API KEY")
    parser.add_argument("--regexes", help="Regex to find environmental samples",
                        type=json.loads )
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--cores', type= int, default=1, help="Cores", required=False)
    parser.add_argument("--log", help="Log file", default= "log.log", required=False)
    
    return parser


##################################################
# FUNCTIONS
##################################################
def fix_host(host, mapping, regexes,
             ask_ncbi = True, ncbi_db = None, ask_two_words = True):
   """
   From a given input host name, try to assign taxa recursively

   ::param host: (str) input host name
   ::param existing_map: (v) existing host mapping
   ::param ask: (bool) If TRUE and more than one TSN is found for the species,
         the user is asked for input. If FALSE NA is returned for multiple matches.
   ::param ask_ever: Wheter to iterate recursively if no results
   ::param bold: (bool) Whether to search by scientific name
   ..param api_key: (str) NCBI API Key
   """
   import re
   from rapidfuzz import process
   from utils import mapping_process_host, mapping_find_patterns
   import logging
   
   my_logger = logging.getLogger(__name__)

   # Trim whitespaces, double and single quotes
   host = mapping_process_host(host, my_logger)

  # If empty
   if not host or host == "nan":
      my_logger.debug("Empty host")
      return host, False, "empty"
  
  # If existing mapping
  # If match in existing mapping with >= 95% similarity
   match, score, choice = process.extractOne(host, set(mapping.keys()))
   if match and score >= 95:
      host_name = mapping[match]
      my_logger.debug(f"{host} in existing map: {mapping[match]}")
      return host, host_name, "mapping"
  
  # If patterns
   pattern, host_name, assigned_method = mapping_find_patterns(
      host, regexes, my_logger)
   if pattern:
      return pattern, host_name, assigned_method
  
   if ask_ncbi and ncbi_db:
      # If no results, search name in NCBI taxonomy db using ete3
      ncbi = NCBITaxa_mod()
      name2taxid = ncbi.get_name_translator([host])
      if name2taxid:
         host_taxid= list(name2taxid.values())[0][0]
         host_name = list(ncbi.get_taxid_translator([host_taxid]).values())[0]
         my_logger.debug(f"Match for {host} -> {host_name}")
         return host, host_name, "ete3"

      # Search for exact match in NCBI taxonomy
      host_name = search_ncbi_name(host)
      if host_name:
         # Search in ete3 to unify nomenclature
         name2taxid = ncbi.get_name_translator([host_name])
         host_taxid= list(name2taxid.values())[0][0]
         host_name = list(ncbi.get_taxid_translator([host_taxid]).values())[0]
         return host, host_name, "ncbi_taxonomy"

      # Search for fuzzy match in NCBI taxonomy
      ## Common names
      host_name = search_fuzzy_ncbi_name(
         ncbi_db = ncbi_db, query=host.capitalize(), sim=0.9, 
         type_name='common')
      if host_name:
         return host, host_name.capitalize(), "fuzzy_common"
      ## Scientific names
      host_name = search_fuzzy_ncbi_name(
         ncbi_db = ncbi_db, query=host.capitalize(), sim=0.9, 
         type_name='sci')
      if host_name:
         return host, host_name.capitalize(), "fuzzy_sci"

      ## Synonyms of scientific names
      host_name = search_fuzzy_ncbi_name(
         ncbi_db = ncbi_db, query=host.capitalize(), sim=0.9, 
         type_name='syn')
      if host_name:
         return host, host_name.capitalize(), "fuzzy_syn"
  
   # If word in parenthesis
   parenthesis = r'\(([^)]+)\)'
   match = re.search(parenthesis, host)
   if match:
      host_name_mod= match.group(1)
      return(fix_host(host_name_mod, mapping,
               regexes = regexes, ask_ncbi = ask_ncbi,
               ncbi_db=ncbi_db, ask_two_words=True))
  
   # Try remove patterns and recursively call
   match = re.search(regexes['patterns'], host)
   if match:
      pattern = match.group()
      host_name_mod = re.sub(pattern, "", host)
      my_logger.debug(f"Recursive call for host = {host_name_mod}")
      return(fix_host(host_name_mod, mapping,
               regexes = regexes, ask_ncbi = ask_ncbi,
               ncbi_db=ncbi_db, ask_two_words=True))

   # Use the first two words
   match = re.search(regexes['two_words'], host)
   if match and ask_two_words:
      host_name_mod = match.group()
      my_logger.debug(f"Recursive call for host = {host_name_mod}")
      return(fix_host(host_name_mod, mapping,
               regexes = regexes, ask_ncbi = ask_ncbi,
               ncbi_db=ncbi_db, ask_two_words=False))
     
   # If none of the previous options have been accomplished
   my_logger.debug(f"{host} -> No result")
   return host, False, "not_assigned"
     
   
def assign_procedure(host_list, mapping, regexes, 
                     ask_ncbi =True, ncbi_db = None, 
                     cores= 10, my_logger=None):
  """
  Approach to assign the host name with its correspodent NCBI taxon UID.
    There are different possible procedures:
      - Previous existence: Assign by existing mapping
      - Instantly assignment: Assign by scientific name
      - Substring assigment: Assign using a substring present in the host name
    ::return data.frame containning the assigned host names (new mapping)
    - Two Columns: Input -> Assigned
  """
  import multiprocessing as mp


   # Split in batches and search
  cmd_batches = [(name, mapping, regexes,
                  ask_ncbi, ncbi_db, True) for name in host_list]
  
  with mp.Pool(processes = cores) as pool:
    try:
        # (Host, host_assigned_name, )
        results = pool.starmap(fix_host, cmd_batches)
    except Exception as e:
        my_logger.error("Something when wrong with fix_host")
        my_logger.error(e)
        pool.terminate()
        raise e
  return results

def search_ncbi_name(host_name):
  from utils import run_cmd

  cmd = f"esearch -db taxonomy -query '{host_name}[All names]' | xtract -pattern ENTREZ_DIRECT -element Count"
  cmd, count, stderr = run_cmd(cmd, split=False, shell=True)
  if count:
   if int(count) == 1:
      cmd = f"esearch -db taxonomy -query '{host_name}[All names]' | esummary | xtract -pattern DocumentSummary -element ScientificName"
      cmd, sci_name, stderr = run_cmd(cmd, split=False, shell=True)
      my_logger.debug(f"Match for {host_name} -> {sci_name.strip()}")
      return sci_name.strip()

def search_fuzzy_ncbi_name(ncbi_db, query, sim=0.9, type_name='common'):
  """Return species name from the NCBI database.

  The results are for the best match for name in the NCBI
  database of taxa names, with a word similarity >= `sim`.

  :param ncbi_db: Species name (does not need to be exact).
  :param query: Species name (does not need to be exact).
  :param 0.9 sim: Min word similarity to report a match (from 0 to 1).
  :param type_name: Which type of name (common, spname, synonym)
  """
  from rapidfuzz import process

  match, score, choice = process.extractOne(query, set(ncbi_db[type_name].values()))
  
  if match and score > int(sim*100):
    ncbi = NCBITaxa_mod()
    taxid = [key for key in ncbi_db[type_name] if ncbi_db[type_name][key] == match][0]
    sci_name = list(ncbi.get_taxid_translator([taxid]).values())[0]
    my_logger.debug(f"Fuzzy match for {query} -> {match} ({sci_name} ; score {score})")

    return sci_name

##################################################
# MAIN
##################################################
if __name__ == "__main__":
   ##################################################
   # FUNCTION ETE3

   from ete3 import NCBITaxa
   class NCBITaxa_mod(NCBITaxa):
      def get_database_names(self):
            """Return taxid, species name and match score from the NCBI database.

            The results are for the best match for name in the NCBI
            database of taxa names, with a word similarity >= `sim`.

            :param name: Species name (does not need to be exact).
            :param 0.9 sim: Min word similarity to report a match (from 0 to 1).
            """
            import sqlite3.dbapi2 as dbapi2
            from utils import merge_dictionary

            _db = dbapi2.connect(self.dbfile)
            # Download scientific and common names
            cmd = (f'SELECT taxid, spname, common FROM species')
            result = _db.execute(cmd)
            dump = result.fetchall()
            db_sci = {taxid: {'spname': spname, 'common': common} for taxid, spname, common in dump}
            
            # Download synonyms
            cmd = (f'SELECT taxid, spname FROM synonym')
            result = _db.execute(cmd)
            dump = result.fetchall()
            db_syn = {taxid: {'synonym': spname} for taxid, spname in dump}

            db_names = merge_dictionary(db_sci, db_syn)
            
            return db_names
   ##################################################

   # Args
   ARGS = get_parser().parse_args() 

   # Logging
   if ARGS.debug:
      logger = My_logger(log_filename = ARGS.log, 
                        logger_name = "create_mappings", level="DEBUG" )
   else:
      logger = My_logger(log_filename = ARGS.log, 
                        logger_name = "create_mappings", level="INFO" )
      
   my_logger = logger.get_logger()


   # Set NBCI API KEY
   os.environ["NCBI_API_KEY"] = f'{ARGS.api_key}'

   # NCBI db
   ncbi = NCBITaxa_mod()
   # NOTE: remove the above commented line. Only to speed up proves
   # ncbi.update_taxonomy_database()
   ncbi_dump = ncbi.get_database_names()
   ncbi_names = {
      'common': {key:entry['common'] for key, entry in ncbi_dump.items()},
      'sci': {key:entry['spname'] for key, entry in ncbi_dump.items()},
      'syn': {key:entry['synonym'] for key, entry in ncbi_dump.items() if 'synonym' in entry.keys()}
      }


   #
   ## Here starts
   #

   df_pls = pd.read_csv(ARGS.input_pls, dtype=str)
   df_pls.loc[:,"BIOSAMPLE_Host"] = df_pls["BIOSAMPLE_Host"].astype("string").str.lower()
   df_pls.loc[:, "BIOSAMPLE_IsolationSource"] = df_pls["BIOSAMPLE_IsolationSource"].astype("string").str.lower()
   my_logger.info(f"Plasmid df: {len(df_pls.index)} entries")

   # Unify missing information labels
   df_pls = unify_missing_info(df_pls, 
                               value_missing="nan", 
                               my_logger=my_logger,
                               use_loop=True)
   host_names = list(df_pls["BIOSAMPLE_Host"])
   iso_names = list(df_pls["BIOSAMPLE_IsolationSource"])

   #
   ## Host inference from BIOSAMPLE Hostname
   #

   # Read Host mapping information: Input -> Assigned
   ## For a given input host name, the assigned scientific name
   df_mapping = pd.read_csv(ARGS.host_mapping)
   my_logger.info(f"Host mapping information ({len(df_mapping.index)} entries) from {ARGS.host_mapping} ...")
   ## Created mapping dictionary with the existing assigned Host names
   host_mapping = create_mapping(df_mapping,
                           input_column="Input_host", 
                           output_column="Final_host")

   # Assign scientific name to BIOSAMPLE_Host names
   my_logger.info(f"""
                  Inferring {len(host_names)} scientific names from BIOSAMPLE Host using host 
                  mapping using {ARGS.cores} cores""")
   host_result = assign_procedure(host_list=host_names, mapping=host_mapping, 
                  regexes=ARGS.regexes, 
                  ask_ncbi = True, ncbi_db=ncbi_names, 
                  cores= ARGS.cores, my_logger=my_logger)
   my_logger.info(f"Host_result len = {len(host_result)}")
   my_logger.info("Assign procedure finished")

   # Discard results from existing mapping
   new_entries = [(input_host, assigned_host, assigned_method
                     ) for input_host, assigned_host, assigned_method in host_result if assigned_method not in ['empty', 'mapping']]
   # Create new df with new mapping results
   if new_entries:
      host_mapping_df = pd.DataFrame(data=new_entries, columns=['Input_host', "Assigned_host", "Assigned_method"], dtype="string")
      host_mapping_df.drop_duplicates(subset=["Input_host", "Assigned_host"], inplace=True)
      host_mapping_df["Eval"] = ""
      host_mapping_df["Final_host"] = ""
      host_mapping_df["Version"] = ARGS.version
      host_mapping_df["Source"] = "Host"
      # Join old mapping with new mapping
      joined_host_mapping_df = pd.concat([df_mapping, host_mapping_df])
      joined_host_mapping_df = joined_host_mapping_df.sort_values(by=[
         "Assigned_host", "Input_host", "Assigned_method"])
      my_logger.info(f"{len(host_mapping_df.index)} new entries to host mapping from BIOSAMPLE_Host: {ARGS.outfile_host_mapping}")
   else:
      joined_host_mapping_df = df_mapping
      my_logger.info(f"0 new entries to host mapping from BIOSAMPLE_Host: {ARGS.outfile_host_mapping}")


   #
   ## Host inference from ISOLATION source
   #

   ## Created mapping dictionary with the existing assigned Host names
   new_host_mapping = create_mapping(joined_host_mapping_df, input_column="Input_host",
                           output_column="Assigned_host")
   mapping = {**host_mapping, **new_host_mapping}

   # Assign scientific name to empty BIOSAMPLE_Host
   my_logger.info(f"Host_res={len(host_result)} ; iso_names={len(iso_names)}")
   to_correct = [name for result_host_mapping, name in zip(host_result, iso_names) if result_host_mapping[1] == False]
   my_logger.info(f"""
                  Inferring {len(to_correct)} scientific names from BIOSAMPLE IsolationSource using host 
                  mapping using {ARGS.cores} cores""")
   iso_result = assign_procedure(host_list=to_correct, mapping=mapping, 
                  regexes=ARGS.regexes, 
                  ask_ncbi = True, ncbi_db=ncbi_names, 
                  cores= ARGS.cores, my_logger=my_logger)
   my_logger.info("Assign procedure finished")

   # Discard results from existing mapping
   new_entries = [(input_host, assigned_host, assigned_method
                     ) for input_host, assigned_host, assigned_method in iso_result if assigned_method not in [
                        'empty', 'mapping']]
   # Create new df with new mapping results
   if new_entries:
      iso_mapping_df = pd.DataFrame(data=new_entries, columns=['Input_host', "Assigned_host", "Assigned_method"], dtype="string")
      iso_mapping_df.drop_duplicates(subset=["Input_host", "Assigned_host"], inplace=True)
      iso_mapping_df["Eval"] = ""
      iso_mapping_df["Final_host"] = ""
      iso_mapping_df["Version"] = ARGS.version
      iso_mapping_df["Source"] = "IsolationSource"
      # Join old mapping with new mapping
      joined_iso_mapping_df = pd.concat([joined_host_mapping_df, iso_mapping_df])
      joined_iso_mapping_df = joined_iso_mapping_df.sort_values(by=[
         "Assigned_host", "Input_host", "Assigned_method"])
      joined_iso_mapping_df.to_csv(ARGS.outfile_host_mapping, index=False, quoting=2)
      joined_iso_mapping_df.to_csv(ARGS.outfile_host_mapping_checked, index=False, quoting=2)
      my_logger.info(f"{len(iso_mapping_df.index)} new entries to host mapping from BIOSAMPLE_ISolationSource: {ARGS.outfile_host_mapping}")
      my_logger.info(f"New mapping file created: {ARGS.outfile_host_mapping}")
   else:
      joined_host_mapping_df.to_csv(ARGS.outfile_host_mapping, index=False, quoting=2)
      joined_host_mapping_df.to_csv(ARGS.outfile_host_mapping_checked, index=False, quoting=2)
      my_logger.info(f"0 new entries host mapping from BIOSAMPLE_ISolationSource: {ARGS.outfile_host_mapping}")

   my_logger.info(
      f"""Please, now inspect manually the mapping file:
      - {ARGS.outfile_host_mapping_checked}
         1) Select the new entries of the mapping (filter by version)
         2) Select the entries different from the Assigned_methods "ncbi_taxonomy" and "ete3" (we assume they are correct)
         3) Evaluate if the Assigned_host is correct:
            - If yes, do nothing
            - If no, write False in the 'Eval' column and write the correct host name in 'Final_host' column
            - If 'Assigned_method' == 'not_assigned':
               - and you think the sample is environmental, write environmetal in the 'Final_host' column
               - If you don't find an appropiate host, write 'unknown' in 'Final_host' column
         4) Please, add new rules to 'fix_host' (process_create_mappings.py) to improve the assignment (if convinient)
         5) Save the document and continue the pipeline :D
            """)
