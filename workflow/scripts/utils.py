#!/usr/bin/python

import re

##################################################
# MISC
##################################################
def split_list(values, size):
    """
    Split given list into chunks of given size
    """
    values = list(values)
    for i in range(0, len(values), size):
        yield values[i:(i + size)]

def split_by_size(input, n):
    """
    Split the input in different groups of size N
    :param input: (int) length of the object
    :param n : (int) size of each group

    :return list of tuples with the start and the end of each group
    """

    ranges = [(i * n, (i + 1) * n)  for i in range((input + n - 1) // n ) ]
    # Adjust the last range to the length of the object
    ranges[-1] = (ranges[-1][0], input)

    return ranges

def asyncio_run_cmd():
    import asyncio
    async def main(args):
        # Create a list of awaitables
        awaitables = [run_cmd2(i) for i in args]
        # gather the awaitables concurrently
        results = await asyncio.gather(*awaitables)
        
        return results


    async def run_cmd2(cmd):

        proc = await asyncio.create_subprocess_shell(
            cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE)

        stdout, stderr = await proc.communicate()

        print(f'[{cmd!r} exited with {proc.returncode}]')
        if stdout:
            print(f'[stdout]\n{stdout.decode()}')
        if stderr:
            print(f'[stderr]\n{stderr.decode()}')
    
    # asyncio.run(main(args=args_cmd))

def run_cmd(cmd, split = True, shell=True):
    import shlex

    """
    Run given CMD
    :param cmd: CMD (str)
    :param split: If the command has to be splitted. Set to False if the command has pipes '|' (bool)
    :param shell: If the command requires shell interpreter. Set to True if the command has pipes '|' (bool)
    :return: cmd (str), status (int), stdout + stderr (str)
    """
    import subprocess

    if split:
        cmd = shlex.split(cmd)

    # check = True, raise an exception if process fail
    try:
        p = subprocess.run(cmd, check=True, 
                        shell=shell, encoding="utf-8",
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE)

        p_stderr = p.stderr
        p_stdout = p.stdout
    
    except subprocess.CalledProcessError as exc:
        print(
            f"Process failed because did not return a successful return code. "
            f"Returned {exc.returncode}\n{exc} \n")
        raise exc

    return cmd, p_stdout, p_stderr

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

def merge_dictionary(dict_1, dict_2):
   dict_3 = {**dict_1, **dict_2}
   for key, value in dict_3.items():
       if key in dict_1 and key in dict_2:
               dict_3[key] = {**value , **dict_1[key]}
   return dict_3

def get_missing_info():
    missing_information = (
        "-","na", "n/a", "n.a.",""," ",
        "missing", "none",  "not applicable",
        "not available", "not collected", "not determined",
        "not recorded",  "unavailable", "unknown",
        "unspecified","unidentified", "null",
        "no", "undocumented", "not defined", "np",
        "not provided", "dead", "alive", "non-human",
        "not reported", "not applied", "dead ?")
    
    return missing_information

def unify_missing_info(df, value_missing="nan", 
                       my_logger=None, use_loop=False):
    missing_information = get_missing_info()
    
    if my_logger:
        my_logger.info(f"Unifying missing information labels... labels for missing information: {missing_information}")
    
    if use_loop:
        for char in missing_information:
            df = df.replace(char, value_missing)
    else:
        df = df.replace(to_replace = missing_information, value=value_missing)
    
    df = df.fillna(value_missing)
    return df

def create_mapping(df, input_column, output_column):
    """
    Create a dictionary from input dataframe
    """
    import numpy as np 
    df = df.loc[:, [input_column, output_column]]
    df = df.set_index(df[input_column]).drop_duplicates()
    mapping = df.to_dict("index")
    mapping = { key:value for key, entry in mapping.items() for colname, value in entry.items()}
    
    if np.nan in mapping.keys():
        del mapping[np.nan]
    
    return mapping

def mapping_process_host(host, my_logger):
    # Trim whitespaces
    host = host.strip().lower()
  
    # Trim possible double quotes
    dd_quotes = r'"'
    host = re.sub(dd_quotes, '', host)
    # Trim possible single quotes
    s_quotes = r"'"
    host = re.sub(s_quotes, '', host)
    
    my_logger.debug(f"Host = {host}")
    return host

def mapping_find_patterns(host, regexes, my_logger):
  import re

  # If environmental sample
  match = re.search(regexes['env'], host)
  if match:
     host_name = match.group()
     my_logger.debug(f"Environmental Match for {host} -> {host_name}")
     return host, host_name, "environmental"
  
  # If Bos taurus synonym
  match = re.search(regexes['bos_taurus'], host)
  if match:
     host_name = "Bos taurus"
     my_logger.debug(f"Match for {host} -> {host_name}")
     return host, host_name, "bos_taurus_synonym"
  
  # If human synonyms
  match = re.search(regexes['human'], host)
  if match:
     host_name = "Homo sapiens"
     my_logger.debug(f"Match for {host} -> {host_name}")
     return host, host_name, "human_synonym"
  
  # If mice synonyms
  match = re.search(regexes['mice'], host)
  if match:
     host_name = "Mus musculus"
     my_logger.debug(f"Match for {host} -> {host_name}")
     return host, host_name, "mouse_synonym"
  
  # If chicken synonyms
  match = re.search(regexes['chicken'], host)
  if match:
     host_name = "Gallus gallus"
     my_logger.debug(f"Match for {host} -> {host_name}")
     return host, host_name, "chicken_synonym"
  
  # If pig synonyms
  match = re.search(regexes['pig'], host)
  if match:
     host_name = "Sus scrofa"
     my_logger.debug(f"Match for {host} -> {host_name}")
     return host, host_name, "pig_synonym"
  
  return None, None, None

def str2timestamp(ts, ts_format='%Y-%m-%d %H:%M:%S'):
    """
    Time stamp string to time stamp object
    """
    import datetime
    return datetime.datetime.strptime(ts, ts_format)

def filter_fasta(ID_list, input_fasta, out_fasta, tmp_file = "tmp_ids.txt",
                 grep_invert = False):
    from os import remove

    with open(tmp_file, 'w') as out:
        out.write("\n".join(str(item) for item in ID_list)) 

    if not grep_invert:
        cmd_fasta = f"seqkit grep -f {tmp_file} {input_fasta} -o {out_fasta}"
    else:
        cmd_fasta = f"seqkit grep -v -f {tmp_file} {input_fasta} -o {out_fasta}"
    
    run_cmd(cmd_fasta,split=True, shell=False)
    cmd_check = f"grep '>' {out_fasta} | wc -l"
    cmd, nrecords, stder = run_cmd(cmd_check, split=False, shell=True)

    remove(tmp_file)
    
    return int(nrecords)

