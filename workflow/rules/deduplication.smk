##################################################
# SEQUENCES DEDUPLICATION
##################################################

rule seqkit_rmdup_seq:
    input:
        fastx=rules.fasta_join.output,
    output:
        fastx=temp(join(OUTDIR, "filtering/deduplication/",
            "seqkit_dedup.fasta")),
        dup_num=join(OUTDIR, "filtering/deduplication/",
            "pls_dup.txt"),
        dup_seqs=temp(join(OUTDIR, "filtering/deduplication/",
            "pls_dupseq.txt")),
    log:
        join(OUTDIR, "filtering/deduplication/", "deduplication.log"),
    benchmark:
        join(OUTDIR, "filtering/deduplication/", "deduplication.bench")
    params:
        command="rmdup",
        extra="--by-seq",
    threads: 4
    wrapper:
        "v3.8.0/bio/seqkit"


rule deduplication:
    input:
        fasta = rules.fasta_join.output,
        fasta_tmp = rules.seqkit_rmdup_seq.output.fastx,
        dup_ids = rules.seqkit_rmdup_seq.output.dup_num,
        dup_seqs_tmp = rules.seqkit_rmdup_seq.output.dup_seqs,
        ## NABT
        nucc = rules.filter_metadata.output.nucc,
        ass = rules.filter_metadata.output.ass,
        bio = rules.filter_metadata.output.bio,
        tax = rules.filter_metadata.output.tax
    output:
        fasta = join(OUTDIR_filtering, "deduplication", "pls_dedup.fasta"),
        nucc_filt = join(OUTDIR_filtering, "deduplication", "nucc_dedup.csv"),
        nucc_identical = join(OUTDIR, "final", "nucc_identical.csv")
    conda: "../envs/seqkit.yaml"
    threads: 1
    log:
        join(OUTDIR_filtering, "deduplication", "deduplication.log"),
    benchmark:
        join(OUTDIR_filtering, "deduplication", "deduplication.bench"),
    script:
        "../scripts/filtering/deduplication.py"

