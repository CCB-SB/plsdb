##################################################
# artifacts
##################################################
rule filter_artifacts:
    input:
        md = rules.chr_filter.output.pls,
        fasta = rules.chr_filter.output.fasta
    output:
        md = join(OUTDIR_filtering, "artifact_filtering.csv"),
        fasta = join(OUTDIR_filtering, "artifact_filtering.fna")
    message:
        "Remove artifacts from plasmid collection"
    params:
        cores = CORES
    log:
        join(OUTDIR_logs, "filtering_artifact.log")
    conda:
        "../../envs/py_env.yml"
    shell:
        """
        python3 scripts/filter_artifacts.py --ifile {input.md} \
            --fasta {input.fasta} --log {log} --out-md {output.md} \
            --out-fasta {output.fasta} --cores {params.cores}
        """