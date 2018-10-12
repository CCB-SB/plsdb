#!/usr/bin/env python3
# This file is from https://bitbucket.org/genomicepidemiology/cge_core_module (1.3.5-2-g08096b1).
# It was slightly modified (see the tag CHANGED).
from __future__ import division
import sys
import os
import time
import random
import re
import subprocess

from Bio.Blast import NCBIXML
from Bio import SeqIO

import collections


class Blaster():
   def __init__(self, inputfile, databases, db_path, out_path='', min_cov=0.6,
                threshold=0.9, blast='blastn', cut_off=True, max_target_seqs=50000):

      min_cov = 100 * float(min_cov)
      threshold = 100 * float(threshold)

      # For alignment data storage
      self.gene_align_query = dict()  # Sequence alignment lines
      self.gene_align_homo = dict()  # Sequence alignment homolog string
      self.gene_align_sbjct = dict()  # Sequence alignment allele string
      self.results = dict()  # Results

      # TODO: Add excluded results to this dictionay
      self.results["excluded"] = dict()

      for db in databases:
         # Adding the path to the database and output
         db_file = "%s/%s" % (db_path, db) # CHANGED by Valentina Galata, 2018.10.11: rm *.fsa
         # tmp_out_path = "%s/tmp" % (out_path)
         # out_file = "%s/out_%s.xml" % (tmp_out_path, db)
         # os.makedirs(tmp_out_path, exist_ok=True)
         # os.chmod(tmp_out_path, 0o775)
         out_file = out_path # CHANGED by Valentina Galata, 2018.10.11: rm *.fsa

         # Running blast
         if os.path.isfile(out_file) and os.access(out_file, os.R_OK):
            print("Found " + out_file + " skipping DB.")
         else:
            cmd = ("%s -db %s -query %s -out %s -outfmt '5'"
                   " -perc_identity  %s -max_target_seqs %s"
                   " -dust 'no'" % (blast, db_file, inputfile,
                       out_file, threshold, max_target_seqs)) # CHANGED by Valentina Galata, 2018.10.11: -db

            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            out, err = process.communicate()

         # Getting the results

         try:
             result_handle = open(out_file, "r")
         except IOError as error:
             sys.exit("Error: BLAST did not run as expected.\n" +
                      "BLAST finished with the following response:" +
                      "\n{}\n{}".format(out.decode("utf-8"),
                                        err.decode("utf-8")))

         blast_records = NCBIXML.parse(result_handle)

         # Declaring variables for saving the results

         # results for each gene
         gene_results = dict()

         # For finding the best hits
         best_hsp = dict()

         # Keeping track of gene split
         gene_split = collections.defaultdict(dict)

         # Making the dicts for sequence outputs
         self.gene_align_query[db] = dict()
         self.gene_align_homo[db] = dict()
         self.gene_align_sbjct[db] = dict()

         # Parsing over the hits and only keeping the best
         for blast_record in blast_records:
            query = blast_record.query
            blast_record.alignments.sort(key=lambda align: (
                -max((len(hsp.query)
                     * (int(hsp.identities) / float(len(hsp.query)))
                     for hsp in align.hsps)))
            )

            for alignment in blast_record.alignments:
               # Setting the e-value as 1 and bit as 0 to get the best
               # HSP fragment
               best_e_value = 1
               best_bit = 0

               for hsp in alignment.hsps:

                  if hsp.expect < best_e_value or hsp.bits > best_bit:
                     best_e_value = hsp.expect
                     best_bit = hsp.bits

                     tmp = alignment.title.split(" ")
                     # sbjct_header = tmp[1]
                     sbjct_header = alignment.title # CHANGED by Valentina Galata, 2018.10.11
                     bit = hsp.bits
                     sbjct_length = alignment.length
                     sbjct_start = hsp.sbjct_start
                     sbjct_end = hsp.sbjct_end
                     gaps = hsp.gaps
                     query_string = str(hsp.query)
                     homo_string = str(hsp.match)
                     sbjct_string = str(hsp.sbjct)
                     contig_name = query.replace(">", "")
                     query_start = hsp.query_start
                     query_end = hsp.query_end
                     HSP_length = len(query_string)
                     perc_ident = int(hsp.identities) / float(HSP_length) * 100
                     strand = 0
                     coverage = ((int(HSP_length) - int(gaps))
                                 / float(sbjct_length))
                     perc_coverage = (((int(HSP_length) - int(gaps))
                                       / float(sbjct_length)) * 100)

                     # cal_score is later used to select the best hit
                     cal_score = perc_ident * coverage

                     hit_id = "%s:%s..%s:%s:%f" % (contig_name, query_start,
                                                   query_end, sbjct_header,
                                                   cal_score)

                     # If the hit is on the other strand
                     if sbjct_start > sbjct_end:
                        tmp = sbjct_start
                        sbjct_start = sbjct_end
                        sbjct_end = tmp

                        query_string = self.reversecomplement(query_string)
                        homo_string = homo_string[::-1]
                        sbjct_string = self.reversecomplement(sbjct_string)
                        strand = 1

                     # Save hit
                     if (cut_off and perc_coverage > 20) or cut_off == False:
                        best_hsp = {'evalue': hsp.expect,
                                    'sbjct_header': sbjct_header,
                                    'bit': bit,
                                    'perc_ident': perc_ident,
                                    'sbjct_length': sbjct_length,
                                    'sbjct_start': sbjct_start,
                                    'sbjct_end': sbjct_end,
                                    'gaps': gaps,
                                    'query_string': query_string,
                                    'homo_string': homo_string,
                                    'sbjct_string': sbjct_string,
                                    'contig_name': contig_name,
                                    'query_start': query_start,
                                    'query_end': query_end,
                                    'HSP_length': HSP_length,
                                    'coverage': coverage,
                                    'cal_score': cal_score,
                                    'hit_id': hit_id,
                                    'strand': strand,
                                    'perc_coverage': perc_coverage
                                    }

               # Saving the result if any
               if best_hsp:
                  save = 1

                  # If there are other gene alignments they are compared
                  if gene_results:
                     tmp_gene_split = gene_split
                     tmp_results = gene_results

                     # Compare the hit results
                     save, gene_split, gene_results = self.compare_results(
                         save, best_hsp, tmp_results, tmp_gene_split)

                  # If the hit is not overlapping with other hit
                  # seqeunces it is kept
                  if save == 1:
                     gene_results[hit_id] = best_hsp

         result_handle.close()

         # If the hit does not cover the entire database reference the
         # missing seqence data are extracted
         keys = list(gene_results.keys())

         for hit_id in keys:
            hit = gene_results[hit_id]

            # Calculate possible split gene coverage
            perc_coverage = hit['perc_coverage']

            if(hit['sbjct_header'] in gene_split
               and len(gene_split[hit['sbjct_header']]) > 1):

               # Calculate new length
               new_length = self.calculate_new_length(gene_split, gene_results,
                                                      hit)
               hit['split_length'] = new_length

               # Calculate new coverage
               perc_coverage = new_length / float(hit['sbjct_length']) * 100

            # If the hit is above the minimum length threshold it is kept
            if perc_coverage >= min_cov:
               if hit['coverage'] == 1:
                  self.gene_align_query[db][hit_id] = hit['query_string']
                  self.gene_align_homo[db][hit_id] = hit['homo_string']
                  self.gene_align_sbjct[db][hit_id] = hit['sbjct_string']
               elif hit['coverage'] != 1:
                  found = False # CHANGED by Valentina Galata, 2018.10.11
                  # Getting the whole database sequence
                  for seq_record in SeqIO.parse(db_file, "fasta"):
                     # if seq_record.description.replace(" ", "") == hit['sbjct_header'].replace(" ", ""):
                     if seq_record.id == hit['sbjct_header'].split(" ")[0]: # CHANGED by Valentina Galata, 2018.10.11
                        found = True # CHANGED by Valentina Galata, 2018.10.11
                        start_seq = str(seq_record.seq)[:int(hit["sbjct_start"])-1]
                        end_seq = str(seq_record.seq)[int(hit["sbjct_end"]):]
                        self.gene_align_sbjct[db][hit_id] = start_seq + hit['sbjct_string'] + end_seq
                        #self.gene_align_sbjct[db][hit_id] = str(seq_record.seq)
                        break
                  assert found, 'Not found: {}'.format(hit['sbjct_header']) # CHANGED by Valentina Galata, 2018.10.11

                  # Getting the whole contig to extract extra query seqeunce
                  contig = ''
                  for seq_record in SeqIO.parse(inputfile, "fasta"):
                     if seq_record.description.replace(" ", "") == hit['contig_name'].replace(" ", ""):
                        contig = str(seq_record.seq)
                        break

                  # Extract extra sequence from query
                  query_seq, homo_seq = self.get_query_align(hit, contig)

                  # Saving the new alignment sequences
                  self.gene_align_query[db][hit_id] = query_seq
                  self.gene_align_homo[db][hit_id] = homo_seq

            else:
               del gene_results[hit_id]
               if hit['sbjct_header'] in gene_split:
                  del gene_split[hit['sbjct_header']]

         # Save the database result
         if gene_results:
            self.results[db] = gene_results
         else:
            self.results[db] = "No hit found"

   @staticmethod
   def reversecomplement(seq):
      # Make reverse complement strand
      trans = str.maketrans("ATGC", "TACG")
      return seq.translate(trans)[::-1]

   @staticmethod
   def compare_results(save, best_hsp, tmp_results, tmp_gene_split):
      """
         Function for comparing hits and saving only the best hit
      """

      # Get data for comparison
      hit_id = best_hsp['hit_id']
      new_start_query = best_hsp['query_start']
      new_end_query = best_hsp['query_end']
      new_start_sbjct = int(best_hsp['sbjct_start'])
      new_end_sbjct = int(best_hsp['sbjct_end'])
      new_score = best_hsp['cal_score']
      new_db_hit = best_hsp['sbjct_header']
      new_contig = best_hsp['contig_name']
      new_HSP = best_hsp['HSP_length']

      # See if the best HSP fragment overlap with another allignment
      # and keep the allignment with the highest score - if the new
      # fragment is not providing new sequence
      keys = list(tmp_results.keys())
      for hit in keys:
         hit_data = tmp_results[hit]
         old_start_query = hit_data['query_start']
         old_end_query = hit_data['query_end']
         old_start_sbjct = int(hit_data['sbjct_start'])
         old_end_sbjct = int(hit_data['sbjct_end'])
         old_score = hit_data['cal_score']
         old_db_hit = hit_data['sbjct_header']
         old_contig = hit_data['contig_name']
         old_HSP = hit_data['HSP_length']

         remove_old = 0

         # If they align to the same gene in the database they are
         # compared
         if new_db_hit == old_db_hit:

            # If the hit provids additional sequence it is kept and
            # the new coverage is saved otherwise the one with the
            # highest score is kept
            if(new_start_sbjct < old_start_sbjct
               or new_end_sbjct > old_end_sbjct):

               # Save the hits as split
               tmp_gene_split[old_db_hit][hit_id] = 1
               if hit not in tmp_gene_split[old_db_hit]:
                  tmp_gene_split[old_db_hit][hit] = 1

            else:
               if new_score > old_score:
                  # Set to remove old hit
                  remove_old = 1

                  # Save a split if the new hit still creats one
                  if(new_db_hit in tmp_gene_split
                     and hit_id not in tmp_gene_split[new_db_hit]):

                     tmp_gene_split[new_db_hit][hit_id] = 1
               else:
                  save = 0

                  # If the old and new hit is not identical the
                  # possible saved gene split for the new hit is
                  # removed
                  if hit_id != hit:
                     if(new_db_hit in tmp_gene_split
                        and hit_id in tmp_gene_split[new_db_hit]):

                        del tmp_gene_split[new_db_hit][hit_id]
                  break

         # If the hits comes form the same part of the contig
         # sequnce but match different genes only the best hit is
         # kept
         if new_contig == old_contig:

            # if the two hits cover the exact same place on the
            # contig only the percentage of identity is compared
            if(old_start_query == new_start_query
               and old_end_query == new_end_query):

               if best_hsp['perc_ident'] > hit_data['perc_ident']:

                  # Set to remove old hit
                  remove_old = 1

                  # Save a split if the new hit still creats one
                  if(new_db_hit in tmp_gene_split
                     and hit_id not in tmp_gene_split[new_db_hit]):

                     tmp_gene_split[new_db_hit][hit_id] = 1

               elif best_hsp['perc_ident'] == hit_data['perc_ident']:
                  # Save both

                  # Save a split if the new hit still creats one
                  if(new_db_hit in tmp_gene_split
                     and hit_id not in tmp_gene_split[new_db_hit]):

                     tmp_gene_split[new_db_hit][hit_id] = 1
               else:
                  save = 0

                  # Remove new gene from gene split if present
                  if(new_db_hit in tmp_gene_split
                     and hit_id in tmp_gene_split[new_db_hit]):

                     del tmp_gene_split[new_db_hit][hit_id]
                  break

            elif((max(old_end_query, new_end_query)
                  - min(old_start_query, new_start_query))
                 <= ((old_end_query - old_start_query)
                     + (new_end_query - new_start_query))):

               if new_score > old_score:

                  # Set to remove old gene
                  remove_old = 1

                  # Save a split if the new hit still creats one
                  if(new_db_hit in tmp_gene_split
                     and hit_id not in tmp_gene_split[new_db_hit]):

                     tmp_gene_split[new_db_hit][hit_id] = 1

               elif new_score == old_score:

                  # If both genes are completly covered the longest
                  # hit is chosen
                  if(int(best_hsp['perc_coverage']) == 100
                     and int(hit_data['perc_coverage']) == 100
                     and new_HSP > old_HSP):

                     # Set to remove old gene
                     remove_old = 1

                  # Save a split if the new hit creats one - both
                  # hits are saved
                  if(new_db_hit in tmp_gene_split
                     and hit_id not in tmp_gene_split[new_db_hit]):

                     tmp_gene_split[new_db_hit][hit_id] = 1
               else:
                  # Remove new gene from gene split if present
                  if(new_db_hit in tmp_gene_split
                     and hit_id in tmp_gene_split[new_db_hit]):

                     del tmp_gene_split[new_db_hit][hit_id]

                  save = 0
                  break

         # Remove old hit if new hit is better
         if remove_old == 1:
            del tmp_results[hit]

            # Remove gene from gene split if present
            if(old_db_hit in tmp_gene_split
               and hit in tmp_gene_split[old_db_hit]):

               del tmp_gene_split[old_db_hit][hit]

      return save, tmp_gene_split, tmp_results

   @staticmethod
   def calculate_new_length(gene_split, gene_results, hit):
      """
         Function for calcualting new length if the gene is split on
         several contigs
      """
      # Looping over splitted hits and calculate new length
      first = 1
      for split in gene_split[hit['sbjct_header']]:

         new_start = int(gene_results[split]['sbjct_start'])
         new_end = int(gene_results[split]['sbjct_end'])

         # Get the first HSP
         if first == 1:
            new_length = int(gene_results[split]['HSP_length'])
            old_start = new_start
            old_end = new_end
            first = 0
            continue
         if new_start < old_start:
            new_length = new_length + (old_start - new_start)
            old_start = new_start

         if new_end > old_end:
            new_length = new_length + (new_end - old_end)
            old_end = new_end

      return(new_length)

   @staticmethod
   def get_query_align(hit, contig):
      """
         Function for extracting extra seqeunce data to the query
         alignment if the full reference length are not covered
      """

      # Getting data needed to extract sequences
      query_seq = hit['query_string']
      homo_seq = hit['homo_string']
      sbjct_start = int(hit['sbjct_start'])
      sbjct_end = int(hit['sbjct_end'])
      query_start = int(hit['query_start'])
      query_end = int(hit['query_end'])
      length = int(hit['sbjct_length'])

      # If the alignment doesn't start at the first position data is
      # added to the begnning
      if sbjct_start != 1:
         missing = sbjct_start - 1

         if(query_start >= missing and hit['strand'] != 1
            or hit['strand'] == 1 and missing <= (len(contig) - query_end)):

            # Getting the query sequence.
            # If the the hit is on the other strand the characters
            # are reversed.
            if hit['strand'] == 1:
               start_pos = query_end
               end_pos = query_end + missing
               chars = contig[start_pos:end_pos]
               chars = Blaster.reversecomplement(chars)
            else:
               start_pos = query_start - missing - 1
               end_pos = query_start - 1
               chars = contig[start_pos:end_pos]

            query_seq = chars + str(query_seq)
         else:
            # Getting the query sequence.
            # If the the hit is on the other strand the characters
            # are reversed.
            if hit['strand'] == 1:
               if query_end == len(contig):
                  query_seq = "-" * missing + str(query_seq)
               else:
                  start_pos = query_end
                  chars = contig[start_pos:]
                  chars = Blaster.reversecomplement(chars)

                  query_seq = ("-" * (missing - len(chars))
                               + chars + str(query_seq))
            elif query_start < 3:
               query_seq = "-" * missing + str(query_seq)
            else:
               end_pos = query_start - 2
               chars = contig[0:end_pos]

               query_seq = ("-" * (missing - len(chars))
                            + chars + str(query_seq))

         # Adding to the homo sequence
         spaces = " " * missing
         homo_seq = str(spaces) + str(homo_seq)

      # If the alignment dosen't end and the last position data is
      # added to the end
      if sbjct_end < length:
         missing = length - sbjct_end

         if(missing <= (len(contig) - query_end) and hit['strand'] != 1
            or hit['strand'] == 1 and query_start >= missing):

            # Getting the query sequence.
            # If the the hit is on the other strand the characters
            # are reversed.
            if hit['strand'] == 1:
               start_pos = query_start - missing - 1
               end_pos = query_start - 1
               chars = contig[start_pos:end_pos]
               chars = Blaster.reversecomplement(chars)
            else:
               start_pos = query_end
               end_pos = query_end + missing
               chars = contig[start_pos:end_pos]

            query_seq = query_seq + chars
         else:
            # If the hit is on the other strand the characters are reversed
            if hit['strand'] == 1:
               if query_start < 3:
                  query_seq = query_seq + "-" * missing
               else:
                  end_pos = query_start - 2
                  chars = contig[0:end_pos]
                  chars = Blaster.reversecomplement(chars)

                  query_seq = (query_seq
                               + chars + "-" * (missing - len(chars)))
            elif query_end == len(contig):
               query_seq = query_seq + "-" * missing
            else:
               start_pos = query_end
               chars = contig[start_pos:]

               query_seq = query_seq + chars + "-" * (missing - len(chars))

         # Adding to the homo sequence
         spaces = " " * int(missing)
         homo_seq = str(homo_seq) + str(spaces)

      return query_seq, homo_seq
