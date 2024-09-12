# NUCCORE
# ##################################################

rule nuccore_queries:
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        esearch_query = config["nuccore"]["esearch"]["query"],
        date_collection_range = config['date_collection_range'],
        batch_size = 8000, # Maximum
        eget_cmd = "| esummary -mode json",
        xtract_cmd = ""
    output:
        DIR = directory(join(OUTDIR, "data_collection/nuccore", 
                f"queries_{config['timestamp']}")),
        batches = expand(
            join(OUTDIR, "data_collection/nuccore", 
                f"queries_{config['timestamp']}", "batch_{batches}.pickle"),
                batches = range(0, config['nuccore']['esearch']['batches']))
    threads: 1
    log:
       join(OUTDIR, "data_collection/nuccore",
            f"queries_{config['timestamp']}/nuccore_queries.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/nuccore/queries"

rule nuccore_api:
    input: 
        pickle = join(rules.nuccore_queries.output.DIR, "batch_{batch}.pickle")
    params:
        ncbi_api = config["eutils"]["api_key"],
        api_file = config["api_key_file"],
        batch_size = 5000
    output:
        pickle = join(OUTDIR, "data_collection/nuccore",
            f"api_{config['timestamp']}", "batch_{batch}.pickle")
    threads: 10 # NOTE: max 10 to prevent too many requests
    benchmark: 
        join(OUTDIR, "data_collection/nuccore",
            f"api_{config['timestamp']}", "batch_{batch}.bench")
    log:
       join(OUTDIR, "data_collection/nuccore",
            f"api_{config['timestamp']}", "batch_{batch}.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/nuccore/api"

rule nuccore_extraction:
    input: 
        pickle = lambda wildcards: expand(rules.nuccore_api.output.pickle, batch=wildcards.batch)
    output:
        csv = join(OUTDIR, "data_collection/nuccore",
                f"extraction_{config['timestamp']}", "batch_{batch}.csv")
    log:
       join(OUTDIR, "data_collection/nuccore",
            f"extraction_{config['timestamp']}", "batch_{batch}.log")
    benchmark:
       join(OUTDIR, "data_collection/nuccore",
            f"extraction_{config['timestamp']}", "batch_{batch}.bench")
    threads: 1
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/nuccore/extraction"

rule nuccore_join:
    input: 
        files = expand(rules.nuccore_extraction.output.csv, 
            batch = range(0, config['nuccore']['esearch']['batches']))
    output: 
        csv = join(OUTDIR, "data_collection/nuccore/",
            f"extraction_{config['timestamp']}", f"nuccore_records.csv")
    threads: 1
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/nuccore/join"

# BIOSAMPLE
##################################################

rule biosample_queries:
    input:
        nuccore_df = rules.nuccore_join.output.csv
    params:
        batch_size = 8000
    output:
        DIR = directory(join(OUTDIR, "data_collection/biosample",
            f"queries_{config['timestamp']}")),
        batches = expand(join(OUTDIR, "data_collection/biosample",
            f"queries_{config['timestamp']}", "batch_{batches}.pickle"),
            batches = range(0, config['biosample']['batches']))
    threads: 1
    benchmark: 
        join(OUTDIR, "data_collection/biosample",
            f"queries_{config['timestamp']}", "biosample_queries.bench")
    log: 
        join(OUTDIR, "data_collection/biosample",
            f"queries_{config['timestamp']}", "biosample_queries.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/biosample/queries"

rule biosample_api:
    input:
        pickle = join(rules.biosample_queries.output.DIR, "batch_{batch}.pickle")
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        batch_size = 500 # if bigger, incorrect parsing
    output:
        pickle = join(OUTDIR, "data_collection/biosample",
            f"api_{config['timestamp']}", "batch_{batch}.pickle"),
        link = join(OUTDIR, "data_collection/biosample",
            f"api_{config['timestamp']}", "batch_link_{batch}.csv"),
    threads: 10 # NOTE: max 10 to prevent too many requests
    benchmark: 
        join(OUTDIR, "data_collection/biosample",
            f"api_{config['timestamp']}", "batch_{batch}.bench")
    log: 
        join(OUTDIR, "data_collection/biosample",
            f"api_{config['timestamp']}", "batch_{batch}.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/biosample/api"

rule biosample_join:
    input:
        pickles = expand(rules.biosample_api.output.pickle,
            batch=range(0, config['biosample']['batches'])),
        links = expand(rules.biosample_api.output.link,
            batch=range(0, config['biosample']['batches'])),
    output:
        pickle = join(OUTDIR, "data_collection/biosample",
            f"api_{config['timestamp']}", "api_results.pickle"),
        links = join(OUTDIR, "data_collection/biosample",
            f"api_{config['timestamp']}", "biosample_linked.csv"),
    run:
    
        import pickle
        import pandas as pd

        stdout_res = []
        for batch in input.pickles:
            sub = pickle.load(open(str(batch), 'rb'))
            stdout_res.extend(sub)
        pickle.dump(stdout_res, open(output.pickle, 'wb'))

        biosample_linked_df = pd.concat([pd.read_csv(str(f)) for f in input.links ], ignore_index=True)
        biosample_linked_df.drop_duplicates(inplace=True)
        biosample_linked_df.to_csv(output.links, index=False)
        

rule biosample_extraction:
    input:
        pickle = rules.biosample_join.output.pickle
    params:
        attributes_plsdb = config['biosample']['attributes'],
    output:
        records = join(OUTDIR, "data_collection/biosample", 
            f"extraction_{config['timestamp']}", "biosample_records.csv"),
        attributes = join(OUTDIR, "data_collection/biosample",
            f"extraction_{config['timestamp']}", "biosample_attributes.csv"),
    threads: 1
    benchmark: 
        join(OUTDIR, "data_collection/biosample",
            f"extraction_{config['timestamp']}", "biosample_extraction.bench")
    log: 
        join(OUTDIR, "data_collection/biosample",
            f"extraction_{config['timestamp']}", "biosample_extraction.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/biosample/extraction"

# ASSEMBLY
##################################################
rule assembly_api:
    input:
        nuccore_df = rules.nuccore_join.output.csv
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        batch_size = 5000
    output:
        pickle = join(OUTDIR, "data_collection/assembly", 
            f"api_{config['timestamp']}", "api_results.pickle"),
        links = join(OUTDIR, "data_collection/assembly", 
            f"api_{config['timestamp']}", "assembly_linked_records.csv")
    threads: 10 # NOTE: max 10 to prevent too many requests
    benchmark: 
        join(OUTDIR, "data_collection/assembly", 
            f"api_{config['timestamp']}", "api_results.bench")
    log:
       join(OUTDIR, "data_collection/assembly",
            f"api_{config['timestamp']}", "api_results.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/assembly/api"
        
rule assembly_extraction:
    input:
        pickle = rules.assembly_api.output.pickle
    output:
        records = join(OUTDIR, "data_collection/assembly", 
            f"extraction_{config['timestamp']}", "assembly_records.csv"),
    threads: 1
    benchmark: 
        join(OUTDIR, "data_collection/assembly", 
            f"extraction_{config['timestamp']}", "assembly_extraction.bench")
    log:
       join(OUTDIR, "data_collection/assembly",
            f"extraction_{config['timestamp']}", "assembly_extraction.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/assembly/extraction"

# TAXONOMY
##################################################
rule taxid:
    input:
        nucc = rules.nuccore_join.output.csv
    output:
        tax = join(OUTDIR, f"data_collection/taxonomy/{config['timestamp']}/",
                "tax_records.csv")
    params:
        taxid_col = "NUCCORE_TaxonID",
        taxid_col2 = "NUCCORE_Description"
    log: join(OUTDIR, f"data_collection/taxonomy/{config['timestamp']}/", "taxonomy.log")
    benchmark: join(OUTDIR, f"data_collection/taxonomy/{config['timestamp']}/", "taxonomy.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/taxonomy"

# NUCCORE-ASSEMBLY-BIOSAMPLE
##################################################

checkpoint join_NABT:
    input:
        nuccore = rules.nuccore_join.output.csv,
        assembly = rules.assembly_extraction.output.records,
        assembly_link = rules.assembly_api.output.links,
        biosample = rules.biosample_extraction.output.records,
        biosample_link = rules.biosample_join.output.links,
        taxonomy = rules.taxid.output.tax,
    output:
        nucc = join(OUTDIR, f"data_collection/database_{config['timestamp']}/NABT",
            "nucc_records.csv"),
        bio = join(OUTDIR, f"data_collection/database_{config['timestamp']}/NABT",
            "bio_records.csv"),
        ass = join(OUTDIR, f"data_collection/database_{config['timestamp']}/NABT",
            "ass_records.csv"),
        tax = join(OUTDIR, f"data_collection/database_{config['timestamp']}/NABT",
            "tax_records.csv"),
        postgres = directory(join(OUTDIR, f"data_collection/database_{config['timestamp']}",
            "plsdb")),
        logfile = "logfile"
    log:
       join(OUTDIR, f"data_collection/database_{config['timestamp']}/NABT",
            "join_NABT.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/utils/join_NABT"


# NUCCORE CHROMOSOMAL DB
# ##################################################

rule nuccorechr_queries:
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        esearch_query = config["nuccore"]["chromosomal"]["esearch"]["query"],
        date_collection_range = config['date_collection_range'],
        batch_size = int(config['nuccore']['chromosomal']['esearch']['batches_size']/10), # Maximum,
        eget_cmd = "| efetch -format docsum",
        xtract_cmd = "| xtract -pattern DocumentSummary -element Id"
    output: 
        DIR = directory(join(OUTDIR, "filtering/chromosomal/nuccore", 
                f"queries_{config['timestamp']}")),
        batches = expand(
            join(OUTDIR, "filtering/chromosomal/nuccore", 
                f"queries_{config['timestamp']}", "batch_{batches}.pickle"),
                batches = range(0, config['nuccore']['chromosomal']['esearch']['batches']))
    log: 
        join(OUTDIR, "filtering/chromosomal/nuccore", 
            f"queries_{config['timestamp']}/nuccoredb_ids.log")
    benchmark:
        join(OUTDIR, "filtering/chromosomal/nuccore", 
            f"queries_{config['timestamp']}/nuccoredb_ids.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/nuccore/queries"

rule nuccorechr_api:
    input:             
        pickle = join(rules.nuccorechr_queries.output.DIR, 
            "batch_{batch}.pickle")
    params:
        ncbi_api = config["eutils"]["api_key"],
        api_file = config["api_key_file"],
        batch_size = config['nuccore']['chromosomal']['esearch']['batches_size']
    output:
        pickle = join(OUTDIR, "filtering/chromosomal/nuccore",
            f"api_{config['timestamp']}", "batch_{batch}.pickle")
    threads: 10 # NOTE: max 10 to prevent too many requests
    benchmark: 
        join(OUTDIR, "filtering/chromosomal/nuccore",
            f"api_{config['timestamp']}", "batch_{batch}.bench")
    log:
       join(OUTDIR, "filtering/chromosomal/nuccore",
            f"api_{config['timestamp']}", "batch_{batch}.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/nuccore/api"
        
rule nuccorechr_extraction:
    input:
        pickle = lambda wildcards: expand(rules.nuccorechr_api.output.pickle,
            batch=wildcards.batch)
    output:
        pickle = join(OUTDIR, "filtering/chromosomal/nuccore", 
                f"extraction_{config['timestamp']}", "batch_{batch}.pickle")
    run:
        import pickle
        
        stdout = pickle.load(open(str(input.pickle), 'rb'))
        res="".join(stdout)
        IDs = [i for i in res.split("\n") if i]
        # print(len(IDs))q
        pickle.dump(IDs, open(str(output.pickle), 'wb'))

rule nuccorechr_fasta_api:
    input:
        pickle = lambda wildcards: expand(rules.nuccorechr_extraction.output.pickle,
                batch=wildcards.batch)
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        batch_size = config['nuccore']['chromosomal']['esearch']['batches_size'],
        format='fasta',
    output: 
        fasta = temp(join(OUTDIR, 
            f"filtering/chromosomal/nuccore/fasta/api_{config['timestamp']}/", 
            "batch_{batch}.fasta"))
    log: 
        join(OUTDIR, 
            f"filtering/chromosomal/nuccore/fasta/api_{config['timestamp']}/", 
            "batch_{batch}.log")
    benchmark:
        join(OUTDIR, 
            f"filtering/chromosomal/nuccore/fasta/api_{config['timestamp']}/", 
            "batch_{batch}.bench")
    threads: 1
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/fasta/api"

rule nuccorechr_fasta_join:
    input: 
        fasta = expand(rules.nuccorechr_fasta_api.output.fasta, 
            batch = range(0, config['nuccore']['chromosomal']['esearch']['batches']))
    output:
        fasta = join(OUTDIR, 
            f"filtering/chromosomal/nuccore/fasta/api_{config['timestamp']}/", 
            "sequences.fasta")
    log: 
        join(OUTDIR, 
            f"filtering/chromosomal/nuccore/fasta/api_{config['timestamp']}/", 
            "fasta_join.log")
    params:
        DIR = lambda wildcards, output: dirname(output.fasta),
    threads: 4
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/fasta/join"