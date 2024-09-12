
rule fasta_queries:
    input: rules.filter_metadata.output.nucc 
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        date_collection_range = config['date_collection_range'],
        batch_size = 8000
    output: 
        DIR = directory(join(OUTDIR,
            f"data_collection/nuccore/fasta/queries_{config['timestamp']}/")),
        batches = expand(
                join(OUTDIR, f"data_collection/nuccore/fasta/queries_{config['timestamp']}/",
                    "batch_{batches}.pickle"),
                batches = range(0, config['nuccore']['fasta']['batches']))
    threads: 1
    log: 
        join(OUTDIR, f"data_collection/nuccore/fasta/queries_{config['timestamp']}/", 
            "fasta_queries.log")
    benchmark:
        join(OUTDIR, f"data_collection/nuccore/fasta/queries_{config['timestamp']}/", 
            "fasta_queries.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/fasta/queries"

rule fasta_api:
    input:
        pickle = join(rules.fasta_queries.output.DIR, "batch_{batch}.pickle")
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        batch_size = 8000,
        format = 'fasta',
    output: 
        fasta = temp(join(OUTDIR, 
            f"data_collection/nuccore/fasta/api_{config['timestamp']}/", 
            "batch_{batch}.fasta"))
    threads: 1
    log:
        join(OUTDIR, 
            f"data_collection/nuccore/fasta/api_{config['timestamp']}/", 
            "batch_{batch}.log")
    benchmark:
        join(OUTDIR, 
            f"data_collection/nuccore/fasta/api_{config['timestamp']}/", 
            "batch_{batch}.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/fasta/api"

rule fasta_join:
    input: 
        fasta = expand(rules.fasta_api.output.fasta, 
            batch = range(0, config['nuccore']['fasta']['batches']))
    output:
        fasta = join(OUTDIR, 
            f"data_collection/nuccore/fasta/api_{config['timestamp']}/", 
            "sequences.fasta")
    log: 
        join(OUTDIR, 
            f"data_collection/nuccore/fasta/api_{config['timestamp']}/", 
            "fasta_join.log")
    params:
        DIR = lambda wildcards, output: dirname(output.fasta),
    threads: 4
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/fasta/join"
