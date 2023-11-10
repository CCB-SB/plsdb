"""
Download rMLST alllele sequences from PubMLST web-site.
NOTE: Requires a user login (!)
WHY:  The access through the RESTful API requires further registration steps and looks laborious.

Uses Python Selenium to login, get the table of files and to download these.
Python Selenium API: https://selenium-python.readthedocs.io/api.html
"""

import re
import sys
import pandas
import logging
import argparse
import os
import shutil
from os.path import join, exists
from time import sleep
from utils_my_logger import My_logger
from utils_decorators import timer, debuger
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By

import logging
logg = logging.getLogger("rmlst_download")

##################################################
# ARGS
##################################################
def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--user', '-u', help="", required=True)
    parser.add_argument('--pw',   '-p', help="", required=True)
    parser.add_argument("--outdir", default = "results", required = False,
                        help="Output directory" )
    parser.add_argument("--outfile", default = "res.fasta", required = False,
                        help="Output file" )
    parser.add_argument("--chrome-dw-dir", default = "downloads", required = False,
                        help="Default Chrome downloads directory" )
    parser.add_argument("--log", default = "log.log", required = False)
    parser.add_argument("--cores", default = 1, type=int, required = False)
    return parser

##################################################
# HELPERS
##################################################
@timer(my_logger=logg)
@debuger(my_logger=logg)
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

##################################################
# HELPERS
##################################################
URL = 'https://pubmlst.org/bigsdb?db=pubmlst_rmlst_seqdef'
TITL = 'Log in'

if __name__ == "__main__":

    # Args
    args = get_arg_parser().parse_args()
    
    # Logger
    logger = My_logger(log_filename = args.log, logger_name = "rmlst_download")
    logger = logger.get_logger()

    # Make sure
    print('Chrome/Firefox will download the files to the default download directory. Make sure it does not cause any problems.')
    answer = input('Go on? [Yes|No]')
    while not re.fullmatch(r'y|yes|n|no', answer, re.IGNORECASE):
        answer = input('Go on? [Yes|No]')
    if not re.fullmatch(r'y|yes', answer, re.IGNORECASE):
        sys.exit(0)
    ####################################
    # Open Chrome
    try:
        driver = webdriver.Chrome()
    except:
        profile = webdriver.FirefoxProfile()
        profile.set_preference("browser.download.folderList", 2)
        profile.set_preference("browser.download.manager.showWhenStarting", False)
        profile.set_preference("browser.download.dir", args.chrome_dw_dir)
        profile.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/x-gzip;text/plain")
        driver=webdriver.Firefox(profile)

    # Go to login page
    driver.get(URL)
    assert driver.title == TITL, 'Unexpected title: {}'.format(driver.title)
    # fill form
    elem_login_user = driver.find_element(By.ID, "user")
    elem_login_pw   = driver.find_element(By.ID, "password_field")
    sleep(3)
    elem_login_user.send_keys(args.user)
    sleep(3)
    elem_login_pw.send_keys(args.pw)
    sleep(3)
    elem_login_pw.send_keys(Keys.RETURN)

    # Go to sequences
    sleep(3)
    driver.find_element(By.ID, 'DOWNLOADS_trigger').click()
    sleep(3)
    link_seqs = driver.find_element(By.LINK_TEXT, 'Allele sequences')
    link_seqs.click()
    # click for getting the list of files
    sleep(3)
    elem_seqs = driver.find_element(By.ID, "j1_2_anchor")
    elem_seqs.click()

    # Get table of files
    #first remove the cookie banner which is in the way
    driver.find_element(By.CLASS_NAME, "cc-compliance").click()
    sleep(3)
    # link_tab = driver.find_element_by_link_text("tab-delimited text") # 2018.03.01: does not work anymore
    link_tab = driver.find_element(By.CSS_SELECTOR, "[title='Download table in tab-delimited text format']")
    link_tab.click()
    sleep(3)
    # download table
    out_table = join(args.outdir, "rmlst_files.tsv")
    cmd = f'wget {driver.current_url} -O {out_table}'
    cmd, cmd_s, cmd_o = run_cmd(cmd)
    assert cmd_s == 0
    # load table
    tab = pandas.read_csv(out_table, sep='\t')
    logger.info('Downloaded table\n{}'.format(tab.head()))

    # Download sequences
    # go back
    sleep(3)
    driver.back()
    elem_seqs = driver.find_element(By.ID, "j1_2_anchor")
    elem_seqs.click()
    sleep(3)
    # download sequences
    index=0
    number_of_fails=0
    dw_files = set()
    for locus in tab["locus"]:
        logger.info('Download {}'.format(locus))
        link_fasta = driver.find_element(By.XPATH, '//a[@href="/bigsdb?db=pubmlst_rmlst_seqdef&page=downloadAlleles&locus={}"]'.format(locus))
        link_fasta.location_once_scrolled_into_view  # scroll down so that the correct link is clicked
        link_fasta.click()
        #check if download successfull
        dw_file = join(args.chrome_dw_dir, f"{locus}.fas")
        dw_files.add(dw_file)         

    # Check if downloand success after
    while(not(all([exists(f) for f in dw_files]))):
       sleep(30)
       complete_dws = sum([exists(f) for f in dw_files])
       logger.info(f'A total of {complete_dws} has been downloaded... Left {len(dw_files) - complete_dws} ')
       sleep(30)
    logger.info('Done. You can close the opened browser.')
    
    # Move files from chrome dw directory to outdir
    logger.info('Moving files to output directory...')
    [shutil.move(f, args.outdir, copy_function = shutil.copytree) for f in dw_files]
    logger.info(f'Done. Files to {args.outdir}.')

    # Concatenate all fasta files in a single one
    logger.info(f"Joinning data into one single fasta ... {args.outfile}")
    os.system(f"cat *.fas >> {args.outfile}")