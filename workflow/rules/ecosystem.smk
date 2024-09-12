# BIOSAMPLE Host
##################################################

rule ecosystem_mapping:
    input:
        bio = rules.filter_metadata.output.bio
    output:
        mapping = join(OUTDIR, "ecosystem", "ecosystem_tags.csv"),
    params:
        api_file = config["api_key_file"],
        ncbi_api = config["eutils"]["api_key"],
        # Mapping params
        version = config['version'],
        previous_version = config['previous_version'],
        previous_mapping = config['ecosystem']['previous_mapping'],
        # Thresholds
        threshold_wratio_ncbi = config['ecosystem']['thresholds']['wratio_ncbi'],
        threshold_tsr_ncbi = config['ecosystem']['thresholds']['tsr_ncbi']
    threads: 1
    log: join(OUTDIR, "ecosystem", "mapping.log")
    benchmark: join(OUTDIR, 'ecosystem', 'mapping.bench')
    wrapper:
        "file:///local/plsdb/master/wrappers/ecosystem"

rule ecosystem_manually:
    input:
        mapping = rules.ecosystem_mapping.output.mapping
    output:
        join(OUTDIR, "ecosystem","manually_inspection.done")
    params:
        mapping_checked = f"../src/ecosystem_tags_{config['version']}_checked.csv"
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/ecosystems_manually.py"

rule ecosystem_infer:
    input:
        bio = rules.filter_metadata.output.bio,
        manual_check = rules.ecosystem_manually.output
    output: join(OUTDIR,"ecosystem","ecosystem_infer.csv")
    params:
        previous_mapping = f"../src/ecosystem_infer_{config['version']}_checked.csv",
        mapping_checked = rules.ecosystem_manually.params.mapping_checked,
    log: join(OUTDIR, "ecosystem","ecosystem_infer.log")
    conda: "../envs/ete4.yaml"
    script:
        "../scripts/ecosystem_infer.py"

rule ecosystem_infer_check:
    input:
        mapping = rules.ecosystem_infer.output[0]
    output:
        join(OUTDIR, "ecosystem","infer_check.done")
    params:
        mapping_checked = f"../src/ecosystem_infer_{config['version']}_checked.csv"
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/ecosystems_infer_check.py"

rule ecosystem_taxid:
    input:
        nucc = rules.ecosystem_infer_check.params.mapping_checked
    output:
        tax = join(OUTDIR, f"ecosystem/taxonomy/", "tax_records.csv")
    params:
        taxid_col = "ECOSYSTEM_taxid",
        taxid_col2 = ""
    log: join(OUTDIR, f"ecosystem/taxonomy/", "taxonomy.log")
    benchmark: join(OUTDIR, f"ecosystem/taxonomy/", "taxonomy.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/ncbi/taxonomy"