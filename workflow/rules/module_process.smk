# rMLST
##################################################

rule process_make_rmlst_blastdb:
    input:
        rules.retrieve_rmlst_data.output.fasta
    output:
        dbs = expand(["{file}.{ext}"], 
            file = join(OUTDIR_retrival, "pubmlst/rmlst", "rmlst.fas"),
            ext=config['blast']['ext'])
    params:
        title = 'rmlst'
    conda:
        "../envs/py_env.yml"
    shell:
        "makeblastdb -in {input} -input_type fasta -dbtype nucl -title '{params.title}' -max_file_sz '3GB' "


rule process_rmlst_blastn:
    input:
        fasta = rules.filter_sequences.output.fasta,
        db = rules.process_make_rmlst_blastdb.input,
        dbs = rules.process_make_rmlst_blastdb.output.dbs
    output:
        join(OUTDIR_process, "rmlst_blastn.rmlst")
    threads: CORES
    params:
        header = config['rmlst']['rmlst_header'],
        ident = config['rmlst']['rmlst_ident'],
        cov = config['rmlst']['rmlst_cov']
    log:
        join(OUTDIR_logs, "process_rmlst_blastn.log")
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python3 scripts/process_rmlst_blastn.py --fasta  {input.fasta} \
        --ofile {output} --db {input.db} --dbs rmlst \
        --cores {threads} --header {params.header} \
        --ident {params.ident} --cov {params.cov} --log {log}
        """


# NUCCORE GC content
##################################################
rule process_calculate_GC:
    input:
        pls = rules.filter_artifacts.output.md,
        fasta = rules.filter_artifacts.output.fasta
    output:
        join(OUTDIR_process, "calculate_gc.csv")
    log:
        join(OUTDIR_logs, "process_calculate_gc.log")
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python scripts/process_calculate_gc.py \
            --input-pls {input.pls} --input-fasta {input.fasta}\
            --outfile {output} --log {log}
        """

# NUCCORE LOCAL DB
##################################################

use rule process_make_rmlst_blastdb as process_make_nuccoredb_blastdb with:
    input:
        rules.retrieve_nuccoredb_seqs.output
    output:
        dbs = expand(["{file}.00.{ext}"], 
            file = join(OUTDIR_retrival, "nuccore_local/",
             "nuccoredb_seqs.fasta"),
            ext=config['blast']['ext'])
    params:
        title = config['rmlst']['nuccore_db_name']


# ABRicate annot
##################################################
# Search abr sequences in fasta
rule process_abricate:
    input:
        fna = rules.filter_artifacts.output.fasta, # Final filtered fasta
        abr = rules.retrieve_abricatedb.output
    output:
        join(OUTDIR_process,"abricate", f"abricate.abr.{{dbs}}") 
    params:
        cov=lambda wildcards: config['abricate']['params'][wildcards.dbs]['cov'],
        ident=lambda wildcards: config['abricate']['params'][wildcards.dbs]['ident'],
        cores = int(CORES/5)
    log:
        join(OUTDIR_logs, f"process_abricate_{{dbs}}.log")
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python3 scripts/process_abricate.py --fna {input.fna} \
            --abr {input.abr} --ofile {output} --cov {params.cov} \
            --ident {params.ident} --env_path $CONDA_PREFIX \
            --wildcard {wildcards.dbs} --log {log} --cores {params.cores}
        """

rule process_join_abricate:
    input:
        expand(["{file}.{rdb}"],
            file = join(OUTDIR_process, "abricate", f"abricate.abr"),
            rdb=config['abricate']['params'].keys())
    output:
        csv = join(OUTDIR_process, f"abricate.csv"),
        tsv = join(OUTDIR_process, f"abricate.tsv")
    params:
        cores=10
    message:
        "Concat ABRicate hits: {input}"
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python scripts/process_join_abricate.py \
            --input {input} \
            --ofile {output.csv} \
            --cores {params.cores}
        """


# BIOSAMPLE Host
##################################################
import json 

rule process_create_host_mapping:
    input:
        pls = rules.process_calculate_GC.output,
        mapping = config['host_mapping']['file']
    output:
        mapping = f"../src/mapping_{config['version']}.csv"
    conda:
        "../envs/ete3.yml"
    log:
        join(OUTDIR_logs, "process_create_host_mapping.log")
    params:
        api_key = config["eutils"]["ncbi_api_key"],
        mapping_checked = f"../src/mapping_{config['version']}_checked.csv",
        cores = 50,
        version = config['version'],
        # Regexs for mapping hosts
        regexes = json.dumps(config['host_mapping']['regexes'])
        
    shell:
        """
        python3 scripts/process_create_host_mapping.py \
        --input-pls {input.pls} \
        --host-mapping {input.mapping} --api-key {params.api_key} \
        --outfile-host-mapping {output.mapping} --version {params.version} \
        --outfile-host-mapping-checked {params.mapping_checked} \
        --regexes '{params.regexes}' \
        --log {log} --cores {params.cores}
        """

rule process_manually_inspect_hosts:
    input:
        mapping = rules.process_create_host_mapping.output.mapping
    output:
        join(OUTDIR_process, "manually_inspection_hosts.done")
    params:
        mapping_checked = f"../src/mapping_{config['version']}_checked.csv"
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python3 scripts/process_manually_inspect_hosts.py \
            --input {params.mapping_checked} \
            --outfile {output}
        """

rule process_infer_host:
    input:
        pls = rules.process_calculate_GC.output,
        manual_check = rules.process_manually_inspect_hosts.output
    output:
        join(OUTDIR_process, "infer_host.csv")
    params:
        mapping_checked = f"../src/mapping_{config['version']}_checked.csv",
        regexes = json.dumps(config['host_mapping']['regexes'])
    log:
        join(OUTDIR_logs, "process_infer_host.log")
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python3 scripts/process_infer_host.py --input-pls {input.pls} \
            --mapping {params.mapping_checked} --log {log} \
            --outfile {output} --regexes '{params.regexes}'
        """


# Ontology
##################################################
rule process_disease_ont:
    input:
        pls = rules.process_infer_host.output,
        dis = rules.retrieve_disease_ont.output
    output:
        join(OUTDIR_process, "disease_ontology.csv")
    params:
        cores=CORES
    conda:
        "../envs/requirements.yml"
    log:
        join(OUTDIR_logs, "process_disease_ont.log")
    shell:
        """
        python3 scripts/process_disease_ont.py \
            -i {input.pls} -d {input.dis} -c {params.cores} \
            -o {output} --log {log}
        """


# Location query
##################################################
rule process_parse_locations:
    input:
        pls = rules.process_disease_ont.output,
        corrections = config['location_mapping']['corrections']
    output:
        pls = join(OUTDIR_process, "parse_locations.csv"),
        mapping = f"../src/locations_{config['version']}.csv"
    params:
        api_keys= config['location_mapping']['google_api'],
        locs=config['location_mapping']['mapping'],
        version = config['version'],
        nohits = f"../src/location_corrections_{config['version']}.csv"
    log:
        join(OUTDIR_logs, "process_parse_locations.log")
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python scripts/process_parse_locations.py \
            --plasmids {input.pls} \
            --ofile {output.pls} --output-mapping {output.mapping} \
            --google-api {params.api_keys} \
            --locs-correction {input.corrections}\
            --output-nohits {params.nohits} \
            --locs {params.locs} --log {log} --version {params.version}
        """

rule process_manually_inspect_locations:
    input:
        rules.process_parse_locations.output.mapping
    output:
        done = join(OUTDIR_process, "manually_inspection_locations.done")
    params:
        corrections = f"../src/location_corrections_{config['version']}.csv"
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python3 scripts/process_manually_inspect_locations.py \
            --input-locations {input} \
            --input-corrections {params.corrections} \
            --outfile {output.done}
        """


# pMLST
#################################################

rule process_make_pmlst_blastdb:
    input:
        rules.retrieve_pmlst_data.output
    output:
        join(OUTDIR, "data_retrival", "pubmlst/pmlst", "pmlst_database.done")
    message:
        "Creating pMLST db"
    conda:
        "../envs/py_env.yml"
    shell:
        "$CONDA_PREFIX/scripts/mlst-make_blast_db && touch {output}"

rule process_pmlst:
    input:
        db = rules.process_make_pmlst_blastdb.output,
        fasta = rules.filter_artifacts.output.fasta, # Final filtered fasta
        pf = rules.process_join_abricate.output.csv
    output:
        join(OUTDIR_process, "pmlst_process.csv")
    params:
        map=config['pmlst']['map'],
        cmd=lambda wildcards:config['pmlst']['cmd'],
        env_path="$CONDA_PREFIX"
    log:
        join(OUTDIR_logs, "process_pmlst.log")
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python3 scripts/process_pmlst.py \
        --db $CONDA_PREFIX/db/blast/mlst.fa --fasta {input.fasta} \
        --pf {input.pf} --ofile {output} \
        --map "{params.map}" --cmd "{params.cmd}" \
        --env_path {params.env_path} --log {log}
        """

# Mash
#################################################

# Create Mash signatures from input
rule process_mash_sketch:
    input:
        rules.filter_artifacts.output.fasta
    output:
        join(OUTDIR_process, "mash","mash_signatures.msh")
    params:
        params=config['mash']['sketch_params']
    conda:
        "../envs/py_env.yml"
    shell:
        """
        mash sketch {params.params} \
            -o $(dirname {output})/$(basename \
            -s .msh {output}) {input}
        """

# Mash dist (all, tab format): compare signatures in input
rule process_mash_dist:
    input:
        rules.process_mash_sketch.output
    output:
        join(OUTDIR_process, "mash","mash_dist.dist")
    params:
        params=config['mash']['dist_params']
    conda:
        "../envs/py_env.yml"
    shell:
        "mash dist {params.params} {input} {input} > {output}"

# Mash dist for highly sim. seq.s (dist cutoff)
rule process_mash_dist_sim:
    input:
        rules.process_mash_sketch.output
    output:
        join(OUTDIR_process, "mash", "mash_dist.distS")
    params:
        params=config['mash']['distS_params']
    conda:
        "../envs/py_env.yml"
    shell:
        "mash dist {params.params} {input} {input} > {output}"

# # Mash dist (only if dist = 0)
# rule mash_dist_zero:
#     input:
#         "{basename}.msh"
#     output:
#         "{basename}.dist0"
#     params:
#         params=config['mash']['dist0_params'],
#         #mash=BIN_MASH
#     message:
#         "Compare signatures in {input}"
#     conda:
#         "../envs/requirements.yml"
#     shell:
#         "mash dist {params.params} {input} {input} > {output}"

# UMAP
#################################################

# Embedding using UMAP on Mash distances
rule process_umap:
    input:
        rules.process_mash_dist.output
    output:
        join(OUTDIR_process, "umap_mash_dist.umap")
    params:
        neighbors=config['umap']['neighbors'],
        components=config['umap']['components'],
        min_dist=config['umap']['min_dist']
    log:
        join(OUTDIR_logs, "process_umap.log")
    conda:
        "../envs/requirements.yml"
    shell:
        """
        python scripts/process_umap.py \
            --input {input} --ofile {output} \
            --log {log} --neighbors {params.neighbors} \
            --components {params.components} \
            --min_dist {params.min_dist}
        """
        

# Info table
##################################################
rule process_infotable:
    input:
        pls = rules.process_parse_locations.output.pls,
        location_checking = rules.process_manually_inspect_locations.output,
        emb = rules.process_umap.output,
        abr = rules.process_join_abricate.output.csv,
        pmlst = rules.process_pmlst.output
    output:
        join(OUTDIR_process, f"infotable_{config['version']}.tsv")
    log:
        join(OUTDIR_logs, "process_infotable.log")
    conda:
        "../envs/py_env.yml"
    shell:
        """
        python scripts/process_infotable.py --pls {input.pls} \
            --emb {input.emb} --abr {input.abr} \
            --pmlst {input.pmlst} --ofile {output} \
            --log {log}
        """
