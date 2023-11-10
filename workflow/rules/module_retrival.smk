
rule retrieve_plasmid_metadata:
    params:
        esearch_query = config["eutils"]["query"]["plasmid"]["esearch_query"],
        ncbi_api = config["eutils"]["ncbi_api_key"],
        email = config["eutils"]["email"],
        batch_size = 5000,
        outdir = join(OUTDIR_retrival, "plasmid"),
        threads = 10 # NOTE: max 10 to prevent too many requests
    output:
        records = expand(join(OUTDIR_retrival,
            "plasmid", "{DB}_records.csv"), 
            DB=["nuccore", "assembly", "biosample", "NAB"]),
        links = expand(join(OUTDIR_retrival,
            "plasmid", "{DB}_linked_records.csv"), 
            DB=["assembly", "biosample"])
    log:
       log = join(OUTDIR_logs, "retrival_metadata.log")
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python3 scripts/retrieve_metadata.py --esearch-query '{params.esearch_query}' \
            --ncbi-api {params.ncbi_api} --email {params.email} \
            --batch-size {params.batch_size} \
            --outdir {params.outdir} --log {log.log} --threads {params.threads}
        """

rule retrieve_plasmid_taxid:
    input:
        merged_csv = rules.retrieve_plasmid_metadata.output.records[3]
    params:
        outdir = join(OUTDIR_retrival, "plasmid")
    output:
        records = expand(join(OUTDIR_retrival, 
            "plasmid", "{DB}_records.csv"), 
            DB=["taxonomy", "NABT"])
    log:
       log = join(OUTDIR_logs, "retrival_taxonomy.log")
    conda:
        "../envs/ete3.yml"
    shell:
        """
        python3 scripts/retrieve_taxid.py --input-file {input.merged_csv}\
            --outdir {params.outdir} --log {log.log}
        """

rule retrieve_fasta:
    input: 
        join(OUTDIR_filtering, "metadata_filtering.csv")
    params:
        ncbi_api = config["eutils"]["ncbi_api_key"],
        email = config["eutils"]["email"],
        batch_size = 8000,
        output_prefix = join(OUTDIR_retrival, "plasmid", "metadata_filtering.fna"),
        threads = 10 # Maximum to prevent too ncbi error (many requests)
    output:
        join(OUTDIR_retrival, "plasmid", "metadata_filtering.fna")
    conda:
        "../envs/py_env.yml"
    log:
        join(OUTDIR_logs, "retrival_fasta_metadata_filtering.log")
    shell:
        """
        python3 scripts/retrieve_fasta.py --input-file {input} \
            --ncbi-api {params.ncbi_api} --email {params.email} \
            --batch-size {params.batch_size}\
            --outfile {params.output_prefix} --log {log} --threads {params.threads}
        """

rule retrieve_pmlst_data:
    params:
        db = config['pmlst']['db'],
        url = config['pmlst']['url'],
        url_schemes = config['pmlst']['url_schemes'],
        outdir = join(os.environ['CONDA_PREFIX'], "/db/pubmlst/"),
        cores = 4 # maximum for this API
    output:
        join(OUTDIR_retrival, "pubmlst/pmlst", "download.done")
    conda:
        "../envs/py_env.yml"
    log: 
        join(OUTDIR_logs, "retrival_pmlst.log")
    shell:
        """
        python3 scripts/retrieve_pmlst.py --db {params.db}\
            --url {params.url} --url-schemes {params.url_schemes}\
            --outdir $CONDA_PREFIX/db/pubmlst/ --outfile {output} \
            --cores {params.cores} --log {log}
        """

rule retrieve_rmlst_data:
    #NOTE: This rule needs to be runned locally (no server) because needs webinterface
    params:
        user = config['rmlst']['user'],
        passwd = config['rmlst']['passwd'],
        chrome_dw_dir = config['rmlst']['default_chrome_dw_dir'],
        outdir = join(OUTDIR_retrival, "pubmlst/rmlst"),
        cores = CORES
    output:
        tsv = join(OUTDIR_retrival, "pubmlst/rmlst", "rmlst_files.tsv"),
        fasta = join(OUTDIR_retrival, "pubmlst/rmlst", "rmlst.fas")
    conda:
        "../envs/py_env.yml"
    log: 
        join(OUTDIR_logs, "retrival_rmlst.log")
    message:
        "Download ribosomal data (allele sequences - FASTA files) from PubMLST "
    shell:
        """
        python3 scripts/retrieve_rmlst.py --user {params.user}\
            --pw {params.passwd} --cores {params.cores} --log {log}\
            --outdir {params.outdir} --chrome-dw-dir {params.chrome_dw_dir}\
            --outfile {output.fasta}
        """

rule retrieve_nuccoredb_ids:
    params:
        query = config['rmlst']['nuccore_db_query'],
        ncbi_api = config["eutils"]["ncbi_api_key"],
        threads = 10
    output:
        join(OUTDIR_retrival, 
        "nuccore_local", "nuccoredb_ids.txt")
    conda:
        "../envs/py_env.yml"
    log:
        join(OUTDIR_logs, "retrival_nuccoredb_ids.log")
    shell:
        """
        python3 scripts/retrieve_nuccoredb_ids.py --ncbi-api {params.ncbi_api}\
            --query '{params.query}' --outfile {output} --log {log}\
            --threads {params.threads}
        """


use rule retrieve_fasta as retrieve_nuccoredb_seqs with:
    input: 
        rules.retrieve_nuccoredb_ids.output
    params:
        ncbi_api = config["eutils"]["ncbi_api_key"],
        email = config["eutils"]["email"],
        batch_size = 8000,
        output_prefix = join(OUTDIR_retrival, "nuccore_local", "nuccoredb_seqs.fasta"),
        threads = 10 # Maximum to prevent too ncbi error (many requests)
    output:
        join(OUTDIR_retrival, "nuccore_local", "nuccoredb_seqs.fasta")
    log:
        join(OUTDIR_logs, "retrival_nuccoredb_seqs.log")

if config['abricate']['replace_getdb']:
    rule retrieve_abricate_getdb:
        input: config['abricate']['getdb_bin']
        output: join(OUTDIR_retrival, "abricate","abricate-get_db.done")
        conda:
            "../envs/py_env.yml"
        shell:
            """
            rsync {input} $CONDA_PREFIX/bin/abricate-get_db && touch {output}
            """
else:
    rule retrieve_abricate_getdb:
        output: join(OUTDIR_retrival,"abricate", "abricate-get_db.done")
        conda:
            "../envs/py_env.yml"
        shell:
            """
            touch {output}
            """

rule retrieve_abricatedb:
    input: rules.retrieve_abricate_getdb.output
    output:
        join(OUTDIR_retrival, "abricate", "download.done")
    params:
        dbs = config['abricate']['dbs']
    message:
        "Update ABRicate databases: {params.dbs}"
    conda:
        "../envs/py_env.yml"
    log:
        join(OUTDIR_logs, "retrival_abricatedb.log")
    shell:
        """
	    python3 scripts/retrieve_abricatedb.py --log {log} \
            --databases {params.dbs} --outfile {output}
        """

rule retrieve_disease_ont:
    output:
        join(OUTDIR_retrival, "disease_ontology", f"ontology.obo")
    params:
        url = config['disease_ontology']['url']
    conda:
        "../envs/py_env.yml"
    log:
        join(OUTDIR_logs, "retrival_disease_ont.log")
    shell:
        """
        (echo "### Downloanding disease ontology... ###\n" &&\
        wget {params.url} -O {output} ) 2> {log}
        """