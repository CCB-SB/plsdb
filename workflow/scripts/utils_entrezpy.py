#! /usr/bin/env/python
# coding: utf-8

## -----
# In-house entrezpy functions
# Author: G. Molano, LA (gonmola@hotmail.es)
# Last modified:
## -----

###########
# IMPORTS #
###########
from datetime import datetime

class MyMultiThreading_ncbi():
    def __init__(self, args_cmd, myfunc, file_tmp = "tmp.batch_0", threads=10, my_logger = None):
        self.args_cmd = args_cmd
        self.myfunc = myfunc
        self.file_tmp = file_tmp
        self.threads = threads
        self.my_logger = my_logger
    
    def run_batches_files(self, ranges_batches):
        """
        Run batches files using MultiThreading approach (I/O functions)
        """
        from concurrent.futures import ThreadPoolExecutor
        from os import remove
        from glob import glob
        from utils import split_by_size
        all_res =[]
        ranges_threads = split_by_size(input=len(ranges_batches), n=self.threads)
        for start, end in ranges_threads:
            with ThreadPoolExecutor(self.threads) as pool:
                try:
                    # issue all tasks to the thread pool
                    futures = [pool.submit(self.myfunc, j[0], j[1], j[2]) for j in self.args_cmd[start:end]]
                    # retrieve all return values in order
                    results = [future.result() for future in futures]
                    all_res.extend(results)
                    
                except Exception as e:
                    self.my_logger.error("Something went wrong with Multithreading")
                    [remove(file) for file in glob(f"{self.file_tmp}*") if self.file_tmp]
                    raise e
            ## Check the stdder
            self.my_logger.info(f"Completed Tasks: {start} - {end}")
            [self.my_logger.debug(f"STDERR: {p_stderr}\n CMD: {cmd}") for cmd, p_stdout, p_stderr in results if p_stderr]
        
        return all_res
        
    def check_batches_files(self, check_bytes=True, check_lines=False,
                            check_fasta=False,
                            batch_nentries = None, last_batch_nentries = None,
                            index_batches = None):
        """
        Check the size of the batches files (in bytes and lines) and rerun
            commands if neccesary
        """
        from concurrent.futures import ThreadPoolExecutor
        from utils import run_cmd, split_by_size
        import os 
        from glob import glob
        
        status_dw = 1
        batches_files = glob(f"{self.file_tmp}*")

        while status_dw != 0:
            fail_items = set()
            for file in batches_files:
                # The order of the batches files is not the same
                # order as obtained by glob func (binary order)
                i = index_batches[file]

                # Number of bytes
                if check_bytes:
                    nbytes = os.path.getsize(file)
                    
                    if nbytes == 0:
                        fail_items.add(i)
                        self.my_logger.error(f"Fail batch file {i} [{nbytes} bytes]")
                
                # Number of Lines
                if check_lines:
                    c, stdout, stderr = run_cmd(f"wc -l {file}", split=True, shell=False)
                    nentries = int(stdout.split(" ")[0])
                if check_fasta:
                    c, stdout, stderr = run_cmd(f"grep '>' {file} | wc -l", split=False, shell=True)
                    nentries = int(stdout)

                if check_lines or check_fasta:
                    if i != len(batches_files)-1:
                      if nentries != batch_nentries:
                          fail_items.add(i)
                          self.my_logger.error(f"Fail batch file {i} [{nentries} entries]. Should be {batch_nentries}")
                    else:
                       if nentries != last_batch_nentries:
                          fail_items.add(i)
                          self.my_logger.error(f"Fail batch file {i} [{nentries} entries]. Should be {last_batch_nentries}")
            if fail_items:
                self.my_logger.info(f"TRY AGAIN: {len(fail_items)} batch files fail")
                self.my_logger.debug(f"TRY AGAIN: Batch files {fail_items}")
                args_fail = [self.args_cmd[i] for i in fail_items]
                # [self.my_logger.debug(f"TRY AGAIN: CMDs {args}") for args in args_fail]

                ranges_threads = split_by_size(input=len(args_fail), n=self.threads)
                for start, end in ranges_threads:
                    with ThreadPoolExecutor(self.threads) as pool:
                        try:
                            # issue all tasks to the thread pool
                            futures = [pool.submit(self.myfunc, j[0], j[1], j[2]) for j in args_fail[start:end]]
                            # retrieve all return values in order
                            results_fail = [future.result() for future in futures]
                        except Exception as e:
                            self.my_logger.error("Something went wrong with Multithreading")
                            [os.remove(file) for file in glob.glob(f"{self.tmp_file}*")]
                            raise e
                    # Check the stdder
                    [self.my_logger.debug(f"STDERR: {p_stderr}\n CMD: {cmd}") for cmd, p_stdout, p_stderr in results_fail if p_stderr]
            else:
                self.my_logger.info("Correct downloading of batches files")
                status_dw = 0
    
    def join_batches_files(self, outfile, header=None, sequentially = True, compress=False):
        from os import remove
        from glob import glob
        from utils import run_cmd
        
        batches_files = glob(f"{self.file_tmp}*")
        if header:
            run_cmd(f"echo '{header}' >> {outfile}", split=False, shell=True) # Add header

        if sequentially:
           self.my_logger.info("Sequential joining of batches files...")
           for file in  batches_files:
              run_cmd(f"cat {file} >> {outfile}", split=False, shell=True)
              remove(file)
        else:
          self.my_logger.info("Non-sequential joining of batches files...")
          run_cmd(f"cat {self.file_tmp}* >> {outfile}", split=False, shell=True)
          [remove(file) for file in batches_files]
        
        self.my_logger.info(f"Done. File created: {outfile}")
        
        if compress:
           self.my_logger.info(f"Compressing {outfile}...")
           run_cmd(f"gzip {outfile}")
           self.my_logger.info(f"Done. File created: {outfile}.gz")


class Docsum_nuccore():
  """Simple data class to store individual sequence Docsum records
  from nuccore db.
  """

  class Subtype:

    def __init__(self, subtype, subname):
      self.strain = None
      self.host = None
      self.country = None
      self.collection = None
      self.collection_date = None

      for i in range(len(subtype)):
        if subtype[i] == 'strain':
          self.stain = subname[i]
        if subtype[i] == 'host':
          self.host = subname[i]
        if subtype[i] == 'country':
          self.country = subname[i]
        if subtype[i] == 'collection_date':
          self.collection_date = subname[i]

  def __init__(self, json_docsum):
    self.error = False
    try:
      self.uid = int(json_docsum['uid'])
      self.term = int(json_docsum['term'])
      self.caption = json_docsum['caption']
      self.title = json_docsum['title']
      self.extra = json_docsum['extra']
      self.gi = int(json_docsum['gi'])
      self.createdate = datetime.strptime(json_docsum['createdate'], '%Y/%m/%d').date()
      self.updatedate = json_docsum['updatedate']
      self.flags = json_docsum['flags']
      self.taxid = int(json_docsum['taxid'])
      self.slen =  int(json_docsum['slen'])
      self.biomol =  json_docsum['biomol']
      self.moltype =  json_docsum['moltype']
      self.topology = json_docsum['topology']
      self.sourcedb = json_docsum['sourcedb']
      self.segsetsize = json_docsum['segsetsize']
      self.projectid = int(json_docsum['projectid'])
      self.genome = json_docsum['genome']
      self.subtype = Docsum_nuccore.Subtype(json_docsum['subtype'].split('|'),
                                    json_docsum['subname'].split('|'))
      self.assemblygi = json_docsum['assemblygi']
      self.assemblyacc = json_docsum['assemblyacc']
      self.tech = json_docsum['tech']
      self.completeness = json_docsum['completeness']
      self.geneticcode = int(json_docsum['geneticcode'])
      self.strand = json_docsum['strand']
      self.organism = json_docsum['organism']
      self.strain = json_docsum['strain']
      self.biosample = json_docsum['biosample']
      self.accessionversion = json_docsum['accessionversion']
    except KeyError:
      print(json_docsum)
      self.error = True

class Docsum_assembly():
  """Simple data class to store individual sequence Docsum records from
  Assembly db
  """
  
  

  class Subtype:

    def __init__(self, subtype, subname):
      self.strain = None
      self.host = None
      self.country = None
      self.collection = None
      self.collection_date = None

      for i in range(len(subtype)):
        if subtype[i] == 'strain':
          self.stain = subname[i]
        if subtype[i] == 'host':
          self.host = subname[i]
        if subtype[i] == 'country':
          self.country = subname[i]
        if subtype[i] == 'collection_date':
          self.collection_date = subname[i]

  def __init__(self, json_docsum):
    

    self.error = False
    try:
      self.uid = int(json_docsum['uid'])
      self.rsuid = json_docsum['rsuid']
      self.gbuid = int(json_docsum['gbuid'])
      self.assemblyaccession = json_docsum['assemblyaccession']
      self.lastmajorreleaseaccession = json_docsum['lastmajorreleaseaccession']
      self.latestaccessiongi = json_docsum['latestaccession']
      self.chainid = int(json_docsum['chainid'])
      self.assemblyname = json_docsum['assemblyname']
      self.ucscname= json_docsum['ucscname']
      self.ensemblname = json_docsum['ensemblname'] 
      self.taxid = int(json_docsum['taxid'])
      self.organism =  json_docsum['organism']
      self.assemblytype = json_docsum['assemblytype']
      self.assemblystatus = json_docsum['assemblystatus']
      self.assemblystatussort = json_docsum['assemblystatussort']
      self.wgs = json_docsum['wgs']
      self.coverage = json_docsum['coverage']
      self.seqreleasedate = datetime.strptime(json_docsum['seqreleasedate'], '%Y/%m/%d %H:%M').date()
      self.submissiondate = datetime.strptime(json_docsum['submissiondate'], '%Y/%m/%d %H:%M').date()
      self.propertylist = json_docsum['propertylist']
      # cross references
      self.biosample = json_docsum['biosampleaccn']
      self.genbankid = json_docsum['synonym']['genbank']
      self.ftppath_genbank = json_docsum['ftppath_genbank']
      self.refseqid = json_docsum['synonym']['refseq']
      self.ftppath_refseq = json_docsum['ftppath_refseq']
      self.ftppath_assembly_rpt = json_docsum['ftppath_assembly_rpt']
    except KeyError:
      print(json_docsum)
      self.error = True


class Docsum_biosample():
  """Simple data class to store individual sequence Docsum records from
  Biosample db
  """

  class Subtype:
    
    def __init__(self, string):
      
      self.geographiclocation = Docsum_biosample.Subtype.find_pattern(
        self,
        pattern = r'geographic location">(.*?)</', string = string)
      self.coordinates = Docsum_biosample.Subtype.find_pattern(
        self,
        pattern = r'latitude and longitude">(.*?)</', string = string)
      self.isolationsource = Docsum_biosample.Subtype.find_pattern(
        self,
        pattern = r'isolation source">(.*?)<', string = string)
      self.host = Docsum_biosample.Subtype.find_pattern(
        self,
        pattern = r'host">(.*?)<', string = string)
      self.hostdisease = Docsum_biosample.Subtype.find_pattern(
        self,
        pattern = r'host disease">(.*?)<', string = string)
      self.collectiondate = Docsum_biosample.Subtype.find_pattern(
        self,
        pattern = r'collection date">(.*?)<', string = string)
      self.sampletype = Docsum_biosample.Subtype.find_pattern(
        self,
        pattern = r'sample type">(.*?)<', string = string)

    def find_pattern(self, pattern, string):
      import re
      match = re.search(pattern, string)
      
      if match:
        return match.group(1)
      else:
        return None


  def __init__(self, json_docsum):
    self.error = False
    try:
      self.uid = int(json_docsum['uid'])
      self.title = json_docsum['title']
      self.accession = json_docsum['accession']
      self.sampledata = Docsum_biosample.Subtype(json_docsum['sampledata'])
    except KeyError:
      print(json_docsum)
      self.error = True
    
