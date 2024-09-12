# Filter viral sequences
#################################################
# NOTE: Virsorter2 presents the higher precision. Benchamarking: https://www.biorxiv.org/content/10.1101/2023.04.26.538077v2.full
# We will be following the SOP/protocol: https://dx.doi.org/10.17504/protocols.io.bwm5pc86

rule virsorter2_db:
    output:
        directory(join(OUTDIR_filtering, "viral", "virsorter2", "database"))
    threads: 30
    benchmark: join(OUTDIR_filtering, "viral", "virsorter2", "database", "virsorter2_db.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/virsorter/db"

rule virsorter2:
    input:
        db = rules.virsorter2_db.output,
        fasta = filtered_fasta
    output: 
        DIR = directory(join(OUTDIR_filtering, "viral", "virsorter2", "intermediate")),
        tsv = join(OUTDIR_filtering, "viral", "virsorter2", "intermediate/final-viral-score.tsv"),
        fasta = join(OUTDIR_filtering, "viral", "virsorter2", "intermediate/final-viral-combined.fa"),
    params:
        extra = ["--keep-original-seq"],
        min_length = config['virsorter2']['min_length'],
        min_score = config['virsorter2']['min_score'],
        groups = config['virsorter2']['groups']
    log: join(OUTDIR_filtering, "viral","virsorter2", "intermediate/virsorter2.log")
    threads: workflow.cores
    benchmark: join(OUTDIR_filtering, "viral", "virsorter2", "intermediate/virsorter2.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/virsorter/run"

# NOTE: Below command gives EORFError (as 31.01.2024).
#       For 2024 Using alternative downloand method but check for future releases
# rule checkv_db:
#     output:
#         directory(join(OUTDIR_filtering, 'viral', 'checkv', 'database'))
#     conda: "../envs/viral.yml"
#     threads: 30
#     benchmark: join(OUTDIR_filtering, "viral", "checkv", "database", "checkv_db.bench")
#     shell:
#         "checkv download_database {output}"

rule checkv_db:
    output:
        DIR = directory(join(OUTDIR_filtering, 'viral', 'checkv', 'database')),
        DIR_db = directory(join(OUTDIR_filtering, 'viral', 'checkv', 'database', 
            f"{config['checkv']['db_version']}"))
    params:
        link = config['checkv']['db_link']
    conda: "../envs/checkv.yml"
    threads: workflow.cores
    benchmark: join(OUTDIR_filtering, "viral", "checkv", "database", "checkv_db.bench")
    shell:
        """
        wget {params.link} -O {output.DIR}.tar.gz && \
            tar -xvf {output.DIR}.tar.gz -C {output.DIR} &&\
            cd {output.DIR_db}/genome_db && \
            diamond makedb --in checkv_reps.faa --db checkv_reps --threads {threads}
        """

rule checkv:
    input:
        db = rules.checkv_db.output.DIR_db,
        fasta = rules.virsorter2.output.fasta
    output: 
        DIR = directory(join(OUTDIR_filtering, "viral", "checkv", "results")),
        contamination = join(OUTDIR_filtering, "viral", "checkv", "results", "contamination.tsv"),
        viruses = join(OUTDIR_filtering, "viral", "checkv", "results", "viruses.fna"),
        proviruses = join(OUTDIR_filtering, "viral", "checkv", "results", "proviruses.fna"),
        combined = join(OUTDIR_filtering, "viral", "checkv", "results", "combined.fna"),
    log: join(OUTDIR_filtering, "viral","checkv", "results/checkv.log")
    threads: workflow.cores
    benchmark: join(OUTDIR_filtering, "viral", "checkv", "results/checkv.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/checkv"

rule virsorter2_2:
    input:
        db = rules.virsorter2_db.output,
        fasta = rules.checkv.output.combined
    output: 
        DIR = directory(join(OUTDIR_filtering, "viral", "virsorter2", "final")),
        fasta = join(OUTDIR_filtering, "viral", "virsorter2", "final/for-dramv/final-viral-combined-for-dramv.fa"),
        tab = join(OUTDIR_filtering, "viral", "virsorter2", "final/for-dramv/viral-affi-contigs-for-dramv.tab"),
    params:
        extra = ["--seqname-suffix-off", "--viral-gene-enrich-off", "--provirus-off", "--prep-for-dramv"],
        min_length = config['virsorter2']['min_length'],
        min_score = config['virsorter2']['min_score'],
        groups = config['virsorter2']['groups']
    log: join(OUTDIR_filtering, "viral","virsorter2", "final/virsorter2.log")
    threads: workflow.cores
    benchmark: join(OUTDIR_filtering, "viral", "virsorter2", "final/virsorter2.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/virsorter/run"

# Manually corrected with: https://github.com/WrightonLabCSU/DRAM/issues/340
rule dramv_db:
    output:
        directory(join(OUTDIR_filtering, 'viral', 'dramv', 'database'))
    params:
        extra = "--skip_uniref"
    threads: workflow.cores
    benchmark: join(OUTDIR_filtering, "viral", "dramv", "database", "dramv_db.bench")
    log: join(OUTDIR_filtering, "viral", "dramv", "database", "dramv_db.log")
    conda: "../envs/dramv.yaml"
    shell:
        """
        DRAM-setup.py prepare_databases --threads {threads} \
            {params.extra} --output_dir {output[0]} 2> {log[0]}
        """

rule dramv_annotate:
    input:
        db = rules.dramv_db.output[0],
        input_fasta = rules.virsorter2_2.output.fasta,
        virsorter_affi_contigs = rules.virsorter2_2.output.tab,
    output:
        DIR = directory(join(OUTDIR_filtering, 'viral', 'dramv', 'annotate'))
    params:
        cmd = "annotate",
        extra = "--keep_tmp_dir --skip_trnascan",
        min_contig_size = config['dramv']['min_contig_size'],
        threads = lambda wildcards, threads: threads
    log: join(OUTDIR_filtering, "viral", "dramv", "annotate", "dramv_annotate.log")
    benchmark: join(OUTDIR_filtering, "viral", "dramv", "annotate", "dramv_annotate.bench")
    threads: workflow.cores
    conda: "../envs/dramv.yaml"
    shell:
        """
        DRAM-v.py {params.cmd} \
            --input_fasta {input.input_fasta} \
            --virsorter_affi_contigs {input.virsorter_affi_contigs} \
            -o {output.DIR}/results \
            {params.extra} --threads {threads} \
            --min_contig_size {params.min_contig_size} 2> {log[0]}
        """

rule dramv_summarize:
    input:
        input_file = join(rules.dramv_annotate.output.DIR, "results/annotations.tsv")
    output:
        DIR = directory(join(OUTDIR_filtering, 'viral', 'dramv', 'summary')),
    params:
        cmd = "distill"
    threads: 1
    log: join(OUTDIR_filtering, "viral", "dramv", "summary", "dramv.log")
    benchmark: join(OUTDIR_filtering, "viral", "dramv", "summary", "dramv_summarize.bench")
    conda: "../envs/dramv.yaml"
    shell:
        """
        DRAM-v.py {params.cmd} \
            --input_file {input.input_file} \
            -o {output.DIR}/results 2> {log[0]}
        """

rule viral_curation:
    input: 
        virsorter2 = rules.virsorter2.output.tsv,
        checkv = rules.checkv.output.contamination,
        dramv = join(rules.dramv_summarize.output.DIR, "results/amg_summary.tsv" )
    output:
        pls = join(OUTDIR_filtering, "viral", "pls_putative.csv"),
        viral_tab = join(OUTDIR, "final", "viral_tab.csv"),
    log: join(OUTDIR_filtering, "viral", "pls_putative.log")
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/filtering/viral_curation.py"