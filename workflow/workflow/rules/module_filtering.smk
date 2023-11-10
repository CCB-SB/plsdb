rule filter_metadata:
    input:
        metadata = rules.retrieve_plasmid_taxid.output.records[1]
    params:
        filter_regex = config['filtering']['dfilter']
    output:
        join(OUTDIR_filtering, "metadata_filtering.csv")
    conda: 
        "../envs/py_env.yml"
    log:
        join(OUTDIR_logs, "filtering_metadata.log")
    shell:
        """
        python3 scripts/filter_metadata.py --input-file {input.metadata} \
            --filter-regex '{params.filter_regex}' --outfile {output} --log {log}
        """

rule filter_sequences:
    input:
        fasta = join(OUTDIR_retrival, "plasmid", "metadata_filtering.fna"),
        metadata = rules.filter_metadata.output
    output:
        fasta = join(OUTDIR_filtering, "sequence_filtering.fna"),
        metadata = join(OUTDIR_filtering, "sequences_filtering.csv")
    conda:
        "../envs/py_env.yml"
    log:
        join(OUTDIR_logs, "filtering_sequences.log")
    shell:
        """
        python3 scripts/filter_sequences.py --fasta-file {input.fasta} \
            --metadata-file {input.metadata} --outfile-seqs {output.fasta}\
            --outfile-md {output.metadata} --log {log}
        """


rule filter_rmlst:
    input:
        pls = rules.filter_sequences.output.metadata,
        fna = rules.filter_sequences.output.fasta,
        rmlst = join(OUTDIR_process, "rmlst_blastn.rmlst"),
        db = expand(["{file}.00.{ext}"], 
            file = join(OUTDIR_retrival, 
            "nuccore_local/", "nuccoredb_seqs.fasta"),
            ext=config['blast']['ext'])
    output:
        md = join(OUTDIR_filtering, "rmlst_filtering.csv"),
        fasta = join(OUTDIR_filtering, "rmlst_filtering.fna")
    conda:
        "../envs/py_env.yml"
    params:
        cores = CORES,
        cutoff = config['rmlst']['rmlst_max_loci'],
        blastn_header = config['rmlst']['blastn_header'],
        blastn_pident = config['rmlst']['blastn_pident'],
        blastn_qcovs = config['rmlst']['blastn_qcovs'],
        blastn_db = join(OUTDIR_retrival, "nuccore_local/",
             "nuccoredb_seqs.fasta")
    log:
        join(OUTDIR_logs, "filtering_rmlst.log")
    shell:
        """
        python3 scripts/filter_rmlst.py --pls {input.pls} \
            --fna {input.fna} --rmlst {input.rmlst} \
            --out-metadata {output.md} --out-seqs {output.fasta} \
            --log {log} --cores {params.cores} \
            --cutoff {params.cutoff} \
            --blastn_header {params.blastn_header} \
            --blastn_pident {params.blastn_pident} \
            --blastn_qcovs {params.blastn_qcovs} \
            --blastn_db {params.blastn_db}  
        """

rule filter_artifacts:
    input:
        md = rules.filter_rmlst.output.md,
        fasta = rules.filter_rmlst.output.fasta
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
        "../envs/py_env.yml"
    shell:
        """
        python3 scripts/filter_artifacts.py --ifile {input.md} \
            --fasta {input.fasta} --log {log} --out-md {output.md} \
            --out-fasta {output.fasta} --cores {params.cores}
        """