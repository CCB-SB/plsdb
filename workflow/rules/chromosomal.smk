
# rmlst candidates
##################################################
rule rmlst_api:
    #NOTE: This rule needs to be runned locally (no server) because needs webinterface
    params:
        api_file = config["api_key_file"],
        api_keyword = config['pubmlst']["rmlst"]["api_key"],
        chrome_dw_dir = config['pubmlst']['rmlst']['default_chrome_dw_dir'],
    output:
        DIR = directory(join(OUTDIR_filtering, "chromosomal/rmlst/db")),
        tsv = join(OUTDIR_filtering, "chromosomal/rmlst/db", "rmlst.tsv"),
        fasta = join(OUTDIR_filtering, "chromosomal/rmlst/db", "rmlst.fas")
    log: join(OUTDIR_filtering, "chromosomal/rmlst/db", "rmlst_api.log")
    benchmark: join(OUTDIR_filtering, "chromosomal/rmlst/db", "rmlst_api.bench")
    threads: 1
    wrapper:
        "file:///local/plsdb/master/wrappers/pubmlst/rmlst"

rule rmlst_blastdb_ref:
    input:
        fasta=rules.rmlst_api.output.fasta
    output:
        multiext(rules.rmlst_api.output.fasta,
            ".ndb", ".nhr", ".nin",
            ".not", ".nsq", ".ntf", ".nto")
    log:
        join(dirname(rules.rmlst_api.output.fasta), "blastdb_ref.log")
    log:
        join(dirname(rules.rmlst_api.output.fasta), "blastdb_ref.bench")
    params:
        "-input_type fasta -blastdb_version 5 -parse_seqids -max_file_sz '3GB' "
    wrapper:
        "v3.8.0/bio/blast/makeblastdb"

# NOTE: rules.rmlst_blastn yields the same results as plsdb:20231103-rules.process_rmlst_blastn
rule rmlst_blastn:
    input:
        query = rules.deduplication.output.fasta,
        blastdb = rules.rmlst_blastdb_ref.output
    output:
        tsv = join(OUTDIR_filtering, "chromosomal/rmlst/blastn", "blastn.tsv")
    log:
       join(OUTDIR_filtering, "chromosomal/rmlst/blastn", "blastn.log")
    benchmark:
       join(OUTDIR_filtering, "chromosomal/rmlst/blastn", "blastn.bench")
    threads: workflow.cores
    params:
        format=f"""6 {' '.join(config['rmlst']['rmlst_header'])}""",
        extra=f"-task blastn -perc_identity {config['rmlst']['rmlst_ident']} -qcov_hsp_perc {config['rmlst']['rmlst_cov']} -evalue {config['rmlst']['rmlst_evalue']} -culling_limit 1"
    wrapper:
        "v3.8.0/bio/blast/blastn"

# nuccore verification
##################################################

# TODO: Put into script and add logging
checkpoint chr_putative_ids:
    input: rules.rmlst_blastn.output.tsv
    output: 
        join(OUTDIR_filtering, "chromosomal/nuccore/", "putative_ids.txt")
    params:
        max_loci = config['rmlst']['rmlst_max_loci']
    run:
        import pandas as pd

        try:
            rmlst = pd.read_table(input[0])
            # Count unique loci per query
            # NOTE: verify if sgi == slocus
            rmlst_hits = rmlst.loc[:, ['qseqid', 'sgi']].value_counts('qseqid')
            IDs = rmlst_hits[rmlst_hits.count > params.max_loci, 'sqeqid']

            with open(output[0], 'w') as out:
                out.write(IDs)
        except Exception as e:
            with open(output[0], 'w') as out:
                out.write("")

# input function for the rule aggregate
def filtered_fasta(wildcards):
    with checkpoints.chr_putative_ids.get().output[0].open() as f:
        if f.read().strip() == "":
            return rules.deduplication.output.fasta
        else:
            return rules.filter_pls_and_fasta.output.fasta

def filtered_pls(wildcards):
    with checkpoints.chr_putative_ids.get().output[0].open() as f:
        if f.read().strip() == "":
            return rules.deduplication.output.nucc_filt
        else:
            return rules.filter_pls_and_fasta.output.csv
# TODO: compare with method plsdb_20231103
# NOTE: not tested because of lack of blast hits in rules.rmlst_blastn.output during 2024
##################################################

rule chr_putative_seqs:
    input:
        fastx=rules.deduplication.output.fasta,
        pattern_file=rules.chr_putative_ids.output[0],
    output:
        fastx=temp(join(OUTDIR_filtering, "chromosomal/nuccore","putative_seqs.fasta")),
    log:
        join(OUTDIR_filtering, "chromosomal/nuccore","putative_seqs.log"),
    params:
        command="grep",
    threads: 4
    wrapper:
        "v3.8.0/bio/seqkit"

rule chr_mash_sketch_ref:
    input: rules.nuccorechr_fasta_join.output.fasta
    output:
        join(OUTDIR_filtering, "chromosomal/nuccore/mash","ref_sketch.msh")
    params:
        cmd = "sketch",
        S = 123, # seed
        k = 21, # k-mer size
        s = 1000, #sketch_size
        o = lambda wildcards, output: splitext(output[0][0]),
        extra = ['-i']
    threads: workflow.cores
    wrapper:
        "file:///local/plsdb/master/wrappers/mash"

rule chr_mash_sketch_query:
    input: rules.chr_putative_seqs.output.fastx
    output:
        join(OUTDIR_filtering, "chromosomal/nuccore/mash","query_sketch.msh")
    params:
        cmd = "sketch",
        S = 123, # seed
        k = 21, # k-mer size
        s = 1000, #sketch_size
        o = lambda wildcards, output: splitext(output[0][0]),
        extra = ['-i']
    threads: workflow.cores
    wrapper:
        "file:///local/plsdb/master/wrappers/mash"

rule chr_mash_dist:
    input:
        ref = rules.chr_mash_sketch_ref.output[0],
        query = rules.chr_mash_sketch_query.output[0],
    output:
        join(OUTDIR_filtering, "chromosomal/nuccore/mash","mash_dist.tsv")
    params:
        d = 0.00123693, # cutoff
        t = lambda wildcards, output: output[0]
    benchmark:
        join(OUTDIR_filtering, "chromosomal/nuccore/mash","mash_dist.bench")
    threads: 20
    wrapper:
        "file:///local/plsdb/master/wrappers/mash"

rule chr_mash_dist_list:
    input:
        tsv = rules.chr_mash_dist.output[0]
    output:
        query_list = join(OUTDIR_filtering, "chromosomal/nuccore/blastn","query_list.txt"),
        ref_list = join(OUTDIR_filtering, "chromosomal/nuccore/blastn","ref_list.txt")
    threads: 1
    shell: 
        """
        cut -f 1 {input.tsv} | uniq > {output.ref_list}
        cut -f 2 {input.tsv} | uniq > {output.query_list}
        """
rule chr_fasta_query:
    input:
        fastx=rules.deduplication.output.fasta,
        pattern_file=rules.chr_mash_dist_list.output.query_list,
    output:
        fastx=temp(join(OUTDIR_filtering, "chromosomal/nuccore/blastn","query.fasta")),
    log:
        join(OUTDIR_filtering, "chromosomal/nuccore/blastn","query.log"),
    params:
        command="grep",
    threads: 4
    wrapper:
        "v3.8.0/bio/seqkit"

rule chr_fasta_ref:
    input:
        fastx=rules.nuccorechr_fasta_join.output.fasta,
        pattern_file=rules.chr_mash_dist_list.output.ref_list,
    output:
        fastx=temp(join(OUTDIR_filtering, "chromosomal/nuccore/blastn","ref.fasta")),
    log:
        join(OUTDIR_filtering, "chromosomal/nuccore/blastn","ref.log"),
    params:
        command="grep",
    threads: 4
    wrapper:
        "v3.8.0/bio/seqkit"

rule chr_blastdb_ref:
    input:
        fasta=rules.chr_fasta_ref.output.fastx
    output:
        multiext(rules.chr_fasta_ref.output.fastx,
            ".ndb", ".nhr", ".nin",
            ".not", ".nsq", ".ntf", ".nto")
    log:
        join(dirname(rules.chr_fasta_ref.output.fastx), "blastdb_ref.log")
    log:
        join(dirname(rules.chr_fasta_ref.output.fastx), "blastdb_ref.bench")
    params:
        "-input_type fasta -blastdb_version 5 -parse_seqids -max_file_sz '3GB' "
    wrapper:
        "v3.8.0/bio/blast/makeblastdb"

rule chr_blastn:
    input:
        query = rules.chr_fasta_query.output.fastx,
        blastdb = rules.chr_blastdb_ref.output[0]
    output:
        tsv = join(OUTDIR_filtering, "chromosomal/nuccore/blastn", "blastn.tsv")
    log:
       join(OUTDIR_filtering, "chromosomal/nuccore/blastn", "blastn.log")
    benchmark:
       join(OUTDIR_filtering, "chromosomal/nuccore/blastn", "blastn.bench")
    threads: workflow.cores
    params:
        format=f"6 '{' '.join(config['rmlst']['blastn_header'])}'",
        extra=f"""-task megablast -perc_identity {config['rmlst']['blastn_pident']} \
            -qcov_hsp_perc {config['rmlst']['blastn_qcovs']} \
            -evalue 0.05 \
            -max_target_seqs 10 -max_hsp 10
            """
    wrapper:
        "v3.8.0/bio/blast/blastn"

# TODO: Put into script and add logging
rule chr_chr_ids:
    input: rules.chr_blastn.output.tsv
    output: 
        join(OUTDIR_filtering, "chromosomal/nuccore/", "chromosomal_ids.txt")
    shell:
        "sed '1d' {input[0]} | cut -f 1  | uniq > {output[0]}"

# TODO: filter pls table and fasta file according to results
rule filter_pls_and_fasta:
    input:
        fasta = rules.deduplication.output.fasta,
        # pls = rules.,
        chr_ids = rules.chr_chr_ids.output[0]
    output:
        fasta=join(OUTDIR_filtering, "chromosomal/plsdb.fasta"),
        csv=join(OUTDIR_filtering, "chromosomal/plsdb.csv"),

    shell:
        "touch {output.fasta}"