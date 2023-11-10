##################################################
# IMPORTS
##################################################

import nltk
nltk.download('stopwords')
from nltk.corpus import stopwords

##################################################
# FUNCTIONS
##################################################

def preprocess_df(df):
    result=df.loc[:,["NUCCORE_ACC","BIOSAMPLE_HostDisease"]]
    header="BIOSAMPLE_HostDisease"
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

def workflow_ont(work_row):
    df2 = preprocess_df(work_row)
    df3 = map_ontologies(df2,ontology_dict={"BIOSAMPLE_HostDisease":DO})
    return df3
##################################################
# ARGS
##################################################

def get_arg_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--disease', '-d', help='Disease ontology (.owl)', required=True)
    parser.add_argument('--ifile', '-i', help="Input file to correct entries", required=True)
    parser.add_argument('--ofile', '-o', help="Output file", required=True)
    parser.add_argument('--cores', '-c', help='Number of cores to use', type=int, default=1)
    parser.add_argument('--log', help="Log File", required=False)
    return parser

##################################################
# MAIN
##################################################
#Note NCBI Taxonomy as ONtology does not provide any good results in contrast to the disease ontology
if __name__ == '__main__':
    import numpy as np
    import pandas as pd
    from parse_ontologies import  Disease_ontology
    import multiprocessing as mp
    from utils_my_logger import My_logger

    # Args
    ARGS = get_arg_parser().parse_args()

    # Logger
    logger = My_logger(log_filename = ARGS.log, logger_name = "process_disease_ont")
    my_logger = logger.get_logger()

    # 
    df = pd.read_csv(ARGS.ifile)
    len_ind = len(df.index)
    my_logger.info("Input plasmid info: {} entries".format(len_ind))
    my_logger.info("Parsing BIOSAMPLE_Host_Disease info from {}".format(ARGS.ifile))
    DO = Disease_ontology(location=ARGS.disease)
    df["BIOSAMPLE_HostDisease_processed"]=np.nan
    my_logger.info("Done")

    #
    my_logger.info("Preprocessing and mapping ontologies in {cores} batches".format(cores=ARGS.cores))
    args = [df[index:index+1] for index in list(range(len_ind))]
    with mp.Pool(processes=ARGS.cores) as pool:
        try:
            result = list(pool.map(workflow_ont, args))
        except Exception as e:
            my_logger.error("Something when wrong with workflow_ont")
            pool.terminate()
            raise e
    my_logger.info("Done")
    my_logger.info("Concatenating the results from mapping ontologies...")
    df3 = pd.concat(result)
    my_logger.info("Concatenated results: {0}".format(len(df3.index)))
    
    my_logger.info("Processing concatenated results...")
    for i in df.index:
            if str(df3["BIOSAMPLE_HostDisease"][i])!="nan":
                if df3["BIOSAMPLE_HostDisease"][i] != "disease" and df3["BIOSAMPLE_HostDisease"][i] != "syndrome":
                    df["BIOSAMPLE_HostDisease_processed"][i]=df3["BIOSAMPLE_HostDisease"][i]

    df2 = pd.read_csv(ARGS.ifile)
    garbage=(
        "-","na", "n/a", "n.a.",""," ",
        "missing", "none",  "not applicable",
        "not available", "not collected", "not determined",
        "not recorded",  "unavailable", "unknown",
        "unspecified","unidentified", "null",
        "no", "undocumented", "not defined", "np",
        "not provided", "dead", "alive", "non-human",
        "not reported", "not applied", "dead ?")
    for j in df2.columns:
        m = df2[j].astype(str).str.lower().isin(garbage)
        df.loc[m, j] = np.NAN
    df2["BIOSAMPLE_HostDisease_processed"]=df["BIOSAMPLE_HostDisease_processed"]
    df2.to_csv(ARGS.ofile, index=False, index_label=False)
    my_logger.info("Done. Created file ({0} entries) {1}".format(len(df2.index), ARGS.ofile))
