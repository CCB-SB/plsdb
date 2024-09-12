use rule chr_mash_sketch_query as imgpr_mash_sketch with:
    input: config['imgpr']['fasta']
    output:
        join(OUTDIR, "comparison/","img_sketch.msh")
    params:
        cmd = "sketch",
        S = 123, # seed
        k = 21, # k-mer size
        s = 1000, #sketch_size
        o = lambda wildcards, output: splitext(output[0][0]),
        extra = ['-i']
    threads: workflow.cores

use rule chr_mash_sketch_query as plsdb_mash_sketch with:
    input: filtered_fasta
    output:
        join(OUTDIR, "comparison/","plsdb_sketch.msh")
    params:
        cmd = "sketch",
        S = 123, # seed
        k = 21, # k-mer size
        s = 1000, #sketch_size
        o = lambda wildcards, output: splitext(output[0][0]),
        extra = ['-i']
    threads: workflow.cores

use rule chr_mash_dist as imgpr_mash_dist with :
    input:
        ref = rules.imgpr_mash_sketch.output[0],
        query = rules.plsdb_mash_sketch.output[0],
    output:
        join(OUTDIR, "comparsion","mash_dist.tsv")
    params:
        d = 0.05, # cutoff
        t = lambda wildcards, output: output[0]
    benchmark:
        join(OUTDIR, "comparsion","mash_dist.bench")
    threads: 20

use rule chr_mash_dist_list as imgpr_mash_dist_list with:
    input:
        tsv = rules.imgpr_mash_dist.output[0]
    output:
        ref_list = join(OUTDIR, "comparison/blastn","imgpr_list.txt"),
        query_list = join(OUTDIR, "comparison/blastn","plsdb_list.txt")

use rule chr_fasta_query as plsdb_fasta_query with:
    input:
        fastx=filtered_fasta,
        pattern_file=rules.imgpr_mash_dist_list.output.query_list,
    output:
        fastx=temp(join(OUTDIR, "comparison/blastn","plsdb.fasta")),
    log:
        join(OUTDIR, "comparison/blastn","plsdb.log")

use rule chr_fasta_ref as imgpr_fasta_ref with:
    input:
        fastx=config['imgpr']['fasta'],
        pattern_file=rules.imgpr_mash_dist_list.output.ref_list,
    output:
        fastx=temp(join(OUTDIR, "comparison/blastn","imgpr.fasta")),
    log:
        join(OUTDIR, "comparison/blastn","imgpr.log")

use rule chr_blastdb_ref as imgpr_blastdb_ref with:
    input:
        fasta=rules.imgpr_fasta_ref.output.fastx
    output:
        multiext(rules.imgpr_fasta_ref.output.fastx,
            ".ndb", ".nhr", ".nin",
            ".not", ".nsq", ".ntf", ".nto")
    log:
        join(dirname(rules.imgpr_fasta_ref.output.fastx), "blastdb_ref.log")
    log:
        join(dirname(rules.imgpr_fasta_ref.output.fastx), "blastdb_ref.bench")

use rule chr_blastn as imgpr_blastn with:
    input:
        query = rules.plsdb_fasta_query.output.fastx,
        blastdb = rules.imgpr_blastdb_ref.output[0]
    output:
        tsv = join(OUTDIR, "comparison/blastn", "blastn.tsv")
    log:
       join(OUTDIR, "comparison/blastn", "blastn.log")
    benchmark:
       join(OUTDIR, "comparison/blastn", "blastn.bench")
    threads: workflow.cores
    params:
        format=f"6 '{' '.join(config['imgpr']['blastn_header'])}'",
        extra=f"""-task megablast -perc_identity {config['imgpr']['blastn_pident']} \
            -qcov_hsp_perc {config['imgpr']['blastn_qcovs']} \
            -evalue 0.05 \
            -max_target_seqs 10 -max_hsp 10
            """