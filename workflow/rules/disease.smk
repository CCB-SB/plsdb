# Ontology
##################################################
rule disease:
    input:
        bio = '../../results/filtering/metadata/bio_filt.csv', # rules.filter_metadata.output.bio
        eco = rules.ecosystem_infer_check.output[0],
    output:
        do_terms = join(OUTDIR, "disease_ontology", "do_terms.csv"),
        symp_terms = join(OUTDIR, "disease_ontology", "symp_terms.csv"),
        mapping = join(OUTDIR, "disease_ontology", "mapping.csv"),
        to_check = join(OUTDIR, "disease_ontology", "to_check.csv"),
    params:
        # Ecosystem mapping containing disease tags
        eco_mapping = rules.ecosystem_infer_check.params.mapping_checked,
        # Mapping params
        version = config['version'],
        previous_version = config['previous_version'],
        previous_mapping = config['disease_ontology']['previous_mapping'],
        # Ontologies
        disease_ont = config['disease_ontology']['DO'],
        symp_ont = config['disease_ontology']['SYMP'],
        # Thresholds
        threshold_wratio_do = config['disease_ontology']['threshold_wratio_do'],
        threshold_wratio_symp = config['disease_ontology']['threshold_wratio_symp'],
        threshold_tsr_do = config['disease_ontology']['threshold_tsr_do'],
        threshold_tsr_symp = config['disease_ontology']['threshold_tsr_symp']
    log:
        join(OUTDIR, "disease_ontology", "disease_ont.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/disease_ont"

rule disease_manually:
    input:
        mapping = rules.disease.output.mapping,
    output:
        join(OUTDIR, "disease_ontology", "manually_inspection.done")
    params:
        mapping_checked = f"../src/disease_tags_{config['version']}_checked.csv"
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/diseases_manually.py"

rule disease_infer:
    input:
        bio = '../../results/filtering/metadata/bio_filt.csv',
        eco = rules.ecosystem_infer_check.output[0],
        manual_check = rules.disease_manually.output
    output: join(OUTDIR,"disease_ontology","disease_infer.csv")
    params:
        eco_mapping = rules.ecosystem_infer_check.params.mapping_checked,
        mapping_checked = rules.disease_manually.params.mapping_checked,
        previous_mapping = "../src/disease_infer_2024_checked.csv",
        # Ontologies
        disease_ont = config['disease_ontology']['DO'],
        symp_ont = config['disease_ontology']['SYMP'],
        # Thresholds
        threshold_wratio_do = config['disease_ontology']['threshold_wratio_do'],
        threshold_wratio_symp = config['disease_ontology']['threshold_wratio_symp'],
        threshold_tsr_do = config['disease_ontology']['threshold_tsr_do'],
        threshold_tsr_symp = config['disease_ontology']['threshold_tsr_symp']
    log: join(OUTDIR, "disease_ontology","disease_infer.log")
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/diseases_infer.py"

rule disease_infer_check:
    input:
        mapping = rules.disease_infer.output[0]
    output:
        join(OUTDIR, "disease_ontology","infer_check.done")
    params:
        mapping_checked = f"../src/disease_infer_{config['version']}_checked.csv"
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/diseases_infer_check.py"