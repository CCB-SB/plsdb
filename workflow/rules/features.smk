
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

##################################################
# NCBI + AMR + ANTISMASH
##################################################

rule features_gbk:
    input:
        fasta = filtered_fasta,
        nucc = filtered_pls,
        bgc = rules.antismash_join.output.tsv,
        genbank = rules.genbank_join.output[0],
        amr = rules.hamronize_dedup.output[0]
    output:
        amr_tab = join(OUTDIR, "final", "amr.tsv"),
        gc_tab = join(OUTDIR, "filtering/metadata/nucc_gc.csv"),
        proteins_tab = join(OUTDIR, "final/proteins.csv"),
        proteins = join(OUTDIR, "final/proteins.fasta"),
        DIR = directory(join(OUTDIR, "final/features/gbk/"))
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/features.py"

rule features_json:
    input: 
        DIR = rules.features_gbk.output.DIR,
        config = "../src/cgview_config.yaml"
    output:
        DIR = directory(join(OUTDIR, "final/features/json")),
        DIR_repo = directory("../scripts/cgview-builder/")
    conda: "../envs/cgview_build.yaml"
    shell:
        """
        mkdir -p {output.DIR}
        git clone https://github.com/stothard-group/cgview-builder.git -b master {output.DIR_repo}
        DIR="{input.DIR}/*"
        for file in $DIR; do
            prefix=$(basename -- "$file" .gbk)
            ruby {output.DIR_repo}/cgview_builder_cli.rb --sequence $file \
                --outfile {output.DIR}/$prefix.json \
                --config {input.config}
        done
        """

# rule cluster_proteins:
#     input:
#         rules.features_gbk.proteins
#     output:
#         join(OUTDIR, "proteins/diamond_clust.tsv")
#     shell:
#         """
#         diamond cluster -d {input[0]} -o {output[0]} \
#             --header --approx-id 
#         """
