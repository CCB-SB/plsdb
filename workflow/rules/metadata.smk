##################################################
# METADATA
##################################################
rule filter_metadata:
    input:
        nucc = rules.join_NABT.output.nucc,
        ass = rules.join_NABT.output.ass,
        bio = rules.join_NABT.output.bio,
        tax = rules.join_NABT.output.tax,
    params:
        filter_regex = config['filter']['description_regex']
    output:
        nucc = join(OUTDIR, "filtering/metadata", "nucc_filt.csv"),
        ass = join(OUTDIR, "final", "assembly.csv"),
        bio = join(OUTDIR, "filtering/metadata", "bio_filt.csv"),
        tax = join(OUTDIR, "filtering/metadata", "tax_filt.csv"),
    conda: "../envs/py_env.yaml"
    log: join(OUTDIR, "filtering/metadata", "filter_metadata.log")
    benchmark: join(OUTDIR, "filtering/metadata", "filter_metadata.bench")
    threads: 1
    script:
        "../scripts/filtering/filter_metadata.py"