# Location query
##################################################
rule locations:
    input:
        bio = rules.filter_metadata.output.bio,
        corrections = config['locations']['corrections']
    output:
        bio_loc = join(OUTDIR, "locations/locations_parsed.csv")
    params:
        mapping = f"../src/locations_{config['version']}.csv",
        api_file = config["api_key_file"],
        google_api = config["locations"]["api_key"],
        locs=config['locations']['mapping'],
        version = config['version'],
        nohits = f"../src/location_corrections_{config['version']}.csv"
    log: join(OUTDIR, "locations/locations_parsed.log")
    benchmark: join(OUTDIR, "locations/locations_parsed.bench")
    threads: 1
    wrapper:
        "file:///local/plsdb/master/wrappers/locations"

rule locations_manually:
    input:
        rules.locations.output
    output:
        done = join(OUTDIR, "locations/locations_manually.done")
    params:
        mapping = rules.locations.params.mapping,
        corrections = f"../src/location_corrections_{config['version']}.csv"
    conda: "../envs/py_env.yaml"
    threads: 1
    script:
        "../scripts/locations_manually.py"

rule locations_infer:
    input:
        loc_check = rules.locations_manually.output.done,
        disease_check = rules.disease_infer_check.output[0],
        bio_loc = rules.locations.output.bio_loc,
    output: join(OUTDIR,"locations","locations_infer.csv")
    params:
        eco_dis_infer = rules.disease_infer_check.params.mapping_checked,
        loc_mapping = rules.locations.params.mapping,
        loc_corrections = rules.locations_manually.params.corrections,
        loc_infer = f"../src/location_infer_{config['version']}.csv"
    log: join(OUTDIR, "locations","locations_infer.log")
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/locations_infer.py"