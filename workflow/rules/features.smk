
##################################################
# NCBI features
##################################################

rule genbank_queries:
    input: filtered_pls
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        date_collection_range = config['date_collection_range'],
        batch_size = 8000
    output: 
        DIR = directory(join(OUTDIR,
            f"data_collection/nuccore/genbank/queries_{config['timestamp']}/")),
        batches = expand(
                join(OUTDIR, f"data_collection/nuccore/genbank/queries_{config['timestamp']}/",
                    "batch_{batches}.pickle"),
                batches = range(0, config['nuccore']['fasta']['batches']))
    threads: 1
    log: 
        join(OUTDIR, f"data_collection/nuccore/genbank/queries_{config['timestamp']}/", 
            "fasta_queries.log")
    benchmark:
        join(OUTDIR, f"data_collection/nuccore/genbank/queries_{config['timestamp']}/", 
            "fasta_queries.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/fasta/queries"

rule genbank_api:
    input:
        pickle = join(rules.genbank_queries.output.DIR, "batch_{batch}.pickle")
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        batch_size = 8000,
        format = 'genbank',
    output: 
        fasta = temp(join(OUTDIR, 
            f"data_collection/nuccore/genbank/api_{config['timestamp']}/", 
            "batch_{batch}.gbk"))
    threads: 1
    log:
        join(OUTDIR, 
            f"data_collection/nuccore/genbank/api_{config['timestamp']}/", 
            "batch_{batch}.log")
    benchmark:
        join(OUTDIR, 
            f"data_collection/nuccore/genbank/api_{config['timestamp']}/", 
            "batch_{batch}.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/fasta/api"

rule genbank_join:
    input: 
        gbk = expand(rules.genbank_api.output.fasta, 
            batch = range(0, config['nuccore']['fasta']['batches']))
    output:
        join(OUTDIR, 
            f"data_collection/nuccore/genbank/api_{config['timestamp']}/", 
            "sequences.gbk")
    shell:
        """
        cat {input.gbk} > {output[0]}
        """

###############################################################################
# Prokaryotic Genome Annotation Pipeline
###############################################################################
rule pgap:
    input: "/local/plsdb/plsdb_2024/workflow/AB011549.2.fasta"
    output: directory(join(OUTDIR, "pgap", "test"))
    container: "docker://ncbi/pgap-utils:2024-07-18.build7555"
    shell:
        """
        ./pgap.py -r -o {output[0]} -g {input[0]} -s 'Escherichia coli'
        """

###############################################################################
# IDENTICAL PROTEIN GROUP
###############################################################################
rule ipg_queries:
    input: join(OUTDIR, "final/proteins.csv")
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        database = 'ipg',
        batch_size = 8000, # Maximum
        eget_cmd = "| esummary -mode json",
        xtract_cmd = ""
    output:
        DIR = directory(join(OUTDIR, "data_collection/ipg", 
                f"queries_{config['timestamp']}")),
        batches = expand(
            join(OUTDIR, "data_collection/ipg", 
                f"queries_{config['timestamp']}", "batch_{batches}.pickle"),
                batches = range(0, config['ipg']['batches']))
    threads: 1
    log:
       join(OUTDIR, "data_collection/ipg",
            f"queries_{config['timestamp']}/ipg_queries.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/ipg/queries"

rule ipg_api:
    input: 
        pickle = join(rules.ipg_queries.output.DIR, "batch_{batch}.pickle")
    params:
        ncbi_api = config["eutils"]["api_key"],
        api_file = config["api_key_file"],
        batch_size = 5000
    output:
        pickle = join(OUTDIR, "data_collection/ipg",
            f"api_{config['timestamp']}", "batch_{batch}.pickle")
    threads: 10 # NOTE: max 10 to prevent too many requests
    benchmark: 
        join(OUTDIR, "data_collection/ipg",
            f"api_{config['timestamp']}", "batch_{batch}.bench")
    log:
       join(OUTDIR, "data_collection/ipg",
            f"api_{config['timestamp']}", "batch_{batch}.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/ipg/api"

rule ipg_extraction:
    input: 
        pickle = lambda wildcards: expand(rules.ipg_api.output.pickle, batch=wildcards.batch)
    output:
        csv = join(OUTDIR, "data_collection/ipg",
                f"extraction_{config['timestamp']}", "batch_{batch}.csv")
    log:
       join(OUTDIR, "data_collection/ipg",
            f"extraction_{config['timestamp']}", "batch_{batch}.log")
    benchmark:
       join(OUTDIR, "data_collection/ipg",
            f"extraction_{config['timestamp']}", "batch_{batch}.bench")
    threads: 1
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/ipg/extraction"

rule ipg_join:
    input: 
        files = expand(rules.ipg_extraction.output.csv, 
            batch = range(0, config['ipg']['batches']))
    output: 
        csv = join(OUTDIR, "data_collection/ipg/",
            f"extraction_{config['timestamp']}", f"ipg_records.csv")
    threads: 1
    run:
        import pickle
        import pandas as pd

        ipg_df = pd.concat([pd.read_csv(str(f)) for f in input.files ], ignore_index=True)
        ipg_df.drop_duplicates(inplace=True)
        ipg_df.to_csv(output.csv, index=False)