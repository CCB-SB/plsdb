rule split_fasta_taxid:
    input:
        fasta = "../../results/filtering/deduplication/pls_dedup.fasta",
        tax = "../../results/filtering/metadata/tax_filt.csv",
        nucc = "../../results/filtering/deduplication/nucc_dedup.csv"
    output:
        bac = join(OUTDIR, "eggnog/input/Bacteria.fa"),
        arc = join(OUTDIR, "eggnog/input/Archaea.fa"),
        DIR = directory(join(OUTDIR, "eggnog/input"))
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/split_fasta_taxid.py"

rule eggnog_mapper_dbs:
    output:
        directory(join(OUTDIR, "eggnog/dbs/hmmer/{name}"))
    params:
        outdir = join(OUTDIR, "eggnog/dbs/"),
        dbname = lambda wildcards: wildcards.name,
        taxid = lambda wildcards: config['eggnog']['hmm_db']['taxid'][wildcards.name],
    threads: 1
    benchmark: join(OUTDIR, "eggnog", "dbs", "download_database_{name}.bench")
    benchmark: join(OUTDIR, "eggnog", "dbs", "download_database_{name}.log")
    wrapper:
        "file:///local/plsdb/master/wrappers/eggnog/db"


rule eggnog_mapper:
    input: 
        fasta = abspath(join(rules.split_fasta_taxid.output.DIR, "{name}.fa")),
        dbs = lambda wildcards: expand(rules.eggnog_mapper_dbs.output, name=wildcards.name)
    output:
        res = abspath(join(OUTDIR, "eggnog/mapper/{name}/eggnog_res.txt"))
    params:
        # Input
        itype = "metagenome",
        qtype = 'seq',
        # Database
        data_dir = lambda wildcards, input: abspath(dirname(dirname(input.dbs[0]))),
        dbtype= 'hmmdb',
        database = lambda wildcards: wildcards.name,
        genepred = "prodigal",
        # Output
        scratch_dir = ".snakemake/shadow/",
        output = lambda wildcards, output: splitext(basename(output.res))[0],
        output_dir = lambda wildcards, output: abspath(dirname(output.res)),
        # Others
        evalue = config['eggnog']['evalue'],
        extra = ['-m hmmer'],
        cpu = lambda wildcards, threads: threads,
        num_servers = 10,
    log: join(OUTDIR, "eggnog/mapper/{name}", "mapper.log")
    benchmark: join(OUTDIR, "eggnog/mapper/{name}", "mapper.bench")
    threads: workflow.cores
    shadow: "shallow"
    wrapper:
        "file:///local/plsdb/master/wrappers/eggnog/mapper"
    
rule eggnog_join:
    input: 
        expand(rules.eggnog_mapper.output.res,
            name = config['eggnog']['hmm_db']['taxid'].keys()
            )
    output:
        join(OUTDIR, "eggnog/mapper/summary/res.txt")
    run:
        # Join
        pass