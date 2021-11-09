import json, sys, csv, pickle

import numpy as np

from parse_ontologies import  Disease_ontology
import pandas as pd
import argparse
import nltk
import tqdm
nltk.download('stopwords')
from nltk.corpus import stopwords
import concurrent.futures
import multiprocessing

def preprocess_df(df):
    result=df[["ACC_NUCCORE","Host_DISEASE"]]
    header="Host_DISEASE"
    for i in df.index:
            val2=str(df[header][i])
            if val2=="nan":
                continue
            val_lower = val2.lower().strip().replace('"', '').replace("'", "")
            clean_string = []
            for word in val_lower.split(" "):
                e = ''.join(e for e in word if e.isalnum() or e != "-")
                if e != "":
                    clean_string.append(e)
            val_split = " ".join(clean_string).strip()
            if val_split == "none" or val_split == "no" or val_split == "na" or val_split == "n/a" or \
                val_split == "null" or val_split == "not" or len(val_split) > 50:
                val2=="nan"
            result[header][i]=val2
    return result

def map_ontologies(df,ontology_dict,treshold=80):
    ontology=list(ontology_dict.keys())[0]
    result=pd.DataFrame(index=df.index, columns=[ontology])
    DO=ontology_dict[ontology]
    for i in df.index:
            val = str(df[ontology][i])
            if val == "nan":
                continue
            querywords = val.split()
            resultwords = [word for word in querywords if word.lower() not in stopwords.words('english')]
            val_split = " ".join(resultwords)
            val_split = val_split.lower().strip()
            fuzzy = DO.get_matching_node(val_split)
            if fuzzy[0] > treshold and fuzzy[1] != "":
                     result[ontology][i]=fuzzy[1]
            elif len(val_split) > 1:
                fuzzy_list = DO.get_matching_node_token_set_ratio(val_split, 95)
                if len(fuzzy_list)==0:
                    continue
                    # get maximum of fuzzy list
                maxim_index=fuzzy_list[0][0]
                for tmp in range(len(fuzzy_list)):
                    if fuzzy_list[tmp][0]>fuzzy_list[maxim_index][0]:
                        maxim_index = tmp
                    # one elem is a triple [[ratio, term, syn], ..]
                    #find best description
                if fuzzy_list[maxim_index][1] != "" and  fuzzy_list[maxim_index][0] >= 95:
                    result[ontology][i]=fuzzy_list[maxim_index][1]
    return result

def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--disease', '-d', help='Disease ontology (.owl)', required=True)
    parser.add_argument('--ifile', '-i', help="Input file to correct entries", required=True)
    parser.add_argument('--ofile', '-o', help="Output file", required=True)
    parser.add_argument('--cores', '-c', help='Number of cores to use', type=int, default=1)
    return parser

#Note NCBI Taxonomy as ONtology does not provide any good results in contrast to the disease ontology
if __name__ == '__main__':
    ARGS = get_arg_parser().parse_args()
    df = pd.read_csv(ARGS.ifile, sep="\t", dtype="unicode")
    len_ind = len(df.index)
    DO = Disease_ontology(location=ARGS.disease)
    df["Host_DISEASE_processed"]=np.nan

    def tmp(index):
        work_row=df[index:index+1]
        df2=preprocess_df(work_row) 
        return map_ontologies(df2,ontology_dict={"Host_DISEASE":DO})
    
    with concurrent.futures.ProcessPoolExecutor(ARGS.cores) as pool:
        result= list(tqdm.tqdm(pool.map(tmp, list(range(len_ind)) ,chunksize=100), total=df.shape[0])) 
    df3=pd.concat(result)
    
    for i in df.index:
            if str(df3["Host_DISEASE"][i])!="nan":
                if df3["Host_DISEASE"][i] != "disease" and df3["Host_DISEASE"][i] != "syndrome":
                    df["Host_DISEASE_processed"][i]=df3["Host_DISEASE"][i]

    df2 = pd.read_csv(ARGS.ifile, sep="\t", dtype="unicode")
    garbage=["-","nan","na", "n/a", "n.a.", "missing", "none",  "not applicable",
             "not available", "not collected", "not determined", "not recorded",
             "unavailable", "unknown", "unspecified","unidentified"]
    for j in df2.columns:
        m = df2[j].str.lower().isin(garbage)
        df.loc[m, j] = np.NAN
    df2["Host_DISEASE_processed"]=df["Host_DISEASE_processed"]
    df2.to_csv(ARGS.ofile, sep='\t', header=True, index=False, index_label=False)
