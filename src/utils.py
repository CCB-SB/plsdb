#!/usr/bin/python

import os
import re
import logging

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
    #if p_status != 0:
        #logging.error("CMD failed: status %d\n%s\n%s" % (p_status, cmd, p_stdout))
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

def download_ncbi_assembly(ftp_path, suffix, odir, file_can_exist=False):
    """
    Download an assembly file from NCBI
    Final URL is created as <ftp_path>/<assembly ID>_<suffix>
    :param ftp_path: FTP directory
    :param suffix: FTP file suffix, e.g. genomic.fna.gz
    :param odir: Where to save the file
    :return Output file if succeeded and None otherwise
    """
    # sample ID from URL
    sample_id = ftp_path.split('/')[-1]
    # output file
    sample_file = os.path.join(odir, "%s_%s" % (sample_id, suffix))
    # complete sample URL
    sample_url = "%s/%s_%s" % (ftp_path, sample_id, suffix)
    # check: output file does not exist
    if not file_can_exist:
        assert not os.path.exists(sample_file), "File exists: %s" % sample_file
    # CMD
    cmd = "curl --output %s --fail %s" % (sample_file, sample_url)
    logging.info('CMD: %s' % cmd)
    # execute, try multiple times
    tried = 0
    while tried < 5:
        cmd, cmd_s, cmd_o = run_cmd(cmd)
        if cmd_s == 0: # OK
            return sample_file, True
        else:
            cmd2, cmd2_s, cmd2_o = run_cmd(cmd + ' --head')
            if re.search("file does not exist", cmd2_o): # nothing to download
                logging.info('NO SUCH FILE: %s' % sample_url)
                rm_file(sample_file)
                return sample_file, False
            tried += 1
    rm_file(sample_file)
    logging.info('FAILED DOWNLOAD: %s' % sample_url)
    return sample_file, False

def is_gzipped(fname):
    """
    Naive way to check whether a file is gzipped
    """
    return re.search('\.gz$', fname) or re.search('\.gzip', fname)

def fasta_records(fname):
    """
    Yield FASTA records from given FASTA file (plain or gzipped)
    """
    from Bio import SeqIO
    if is_gzipped(fname):
        import gzip
        with gzip.open(fname, mode='rt') as ifile:
            for record in SeqIO.parse(ifile, 'fasta'):
                yield record
    else:
        with open(fname, mode='r') as ifile:
            for record in SeqIO.parse(ifile, 'fasta'):
                yield record

def records2fasta(records, fname):
    """
    Save given FASTA records to given file
    """
    import gzip
    from Bio import SeqIO
    if is_gzipped(fname):
        with gzip.open(fname, 'wt') as ofile:
            for record in records:
                SeqIO.write(record, ofile, 'fasta')
    else:
        with open(fname, 'w') as ofile:
            for record in records:
                SeqIO.write(record, ofile, 'fasta')

def seq_is_plasmid(seq_header):
    """
    Whether the given sequence header indicates that the
    corresponding sequnece is a plasmid
    """
    return re.search('plasmid', seq_header, re.IGNORECASE)

def adjust_ncbi_tab_header(df_cols):
    df_cols = [re.sub("#\s+", "", h) for h in df_cols]
    return df_cols

def read_and_concat_tabs(tab_files, check_dups=True, sep='\t', header=0):
    """
    Read in and concat multiple tables/data.frames
    """
    import pandas
    # read and concat
    df = pandas.concat([
        pandas.read_csv(filepath_or_buffer=f, sep=sep, header=header)
        for f in tab_files
    ])
    # check for duplicate rows
    if check_dups:
        dupl = df.duplicated()
        assert sum(dupl) == 0, 'Found %d duplicate rows in %s' % (
            sum(dupl),
            ', '.join(tab_files)
        )
    return df

def read_ncbi_tables(tab_files):
    """
    Read in NCBI tables
    If multiple - concat
    Adjust header
    """
    df = read_and_concat_tabs(tab_files)
    df.columns = adjust_ncbi_tab_header(df.columns)
    return df

def chunk_list(ilist, num=1):
    """
    Split list into given number of chunks
    """
    import math
    ilist = list(ilist)
    n = math.ceil(len(ilist) / num)
    for i in range(0, len(ilist), n):
        yield ilist[i:i + n]
