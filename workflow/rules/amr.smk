OUTDIR_amr = join(OUTDIR, "AMR")
PARTS = config["amr"]['batches']

rule seqkit_split_amr:
    input:
        fastx=filtered_fasta
    output:
        fastx=temp(expand(join(OUTDIR_amr, "split/pls_dedup.{part}.fasta"),
            part=[f"part_{str(i).rjust(3, "0")}" for i in range(1, PARTS+1)]))
    log:
        join(OUTDIR_amr, "split/seqkit_split.log"),
    benchmark:
        join(OUTDIR_amr, "split/seqkit_split.bench"),
    params:
        command="split2",
        extra=f"""--by-part {PARTS} --force \
            --out-dir {join(OUTDIR_amr, "split/")} 
            """
    threads: 4
    conda: "../envs/seqkit.yaml"
    shell:
        """
        seqkit {params.command} {input.fastx} --threads {threads} {params.extra} 2> {log}
        """

# NCBI AMRFinder
#################################################
rule AMRFinderPlus_getdb:
    output:
        directory(join(OUTDIR_amr, "amrfinderplus", "database"))
    benchmark: join(OUTDIR_amr, "amrfinderplus/database/benchmark.txt")
    log: join(OUTDIR_amr, "amrfinderplus/database/log.txt")
    threads: workflow.cores
    wrapper:
        "file:///local/plsdb/master/wrappers/amrfinderplus/db"
    

rule AMRFinderPlus_run:
    input:
        nucleotide = join(dirname(rules.seqkit_split_amr.output.fastx[0]), "pls_dedup.{part}.fasta"),
        database = join(rules.AMRFinderPlus_getdb.output[0], "latest")
    output:
        output = join(OUTDIR_amr, "amrfinderplus/{part}/", "results.txt"),
        nucleotide_output = join(OUTDIR_amr, "amrfinderplus/{part}/", "results_nucleotides.txt")
    params:
        ident_min = (config['amrfinderplus']['ident_min'])/100,
        coverage_min = (config['amrfinderplus']['coverage_min'])/100,
        extra = ['--report_all_equal', '--plus']
    benchmark: join(OUTDIR_amr, "amrfinderplus/{part}/amrfinderplus.bench")
    log: join(OUTDIR_amr, "amrfinderplus/{part}/amrfinderplus.log")
    threads: 1  # Do not scale well with threads
    wrapper:
        "file:///local/plsdb/master/wrappers/amrfinderplus/run"

rule AMRFinderPlus_join:
    input:
        expand(rules.AMRFinderPlus_run.output.output,
            part=[f"part_{str(i).rjust(3, "0")}" for i in range(1, PARTS+1)])
    output:
        tsv = join(OUTDIR_amr, "amrfinderplus/summary", "plsdb.tsv")
    threads: 1
    run:
        import pandas as pd
        df = pd.concat([pd.read_table(str(f)) for f in input ], ignore_index=True)
        df.loc[:, 'Protein identifier'] = 'NA'
        df.to_csv(output.tsv, sep="\t", index=False)

# CARD
#################################################
rule rgi_getdb:
    output: directory(join(OUTDIR_amr, "rgi/db"))
    wrapper:
        "file:///local/plsdb/master/wrappers/rgi/db"

# For faster searching, use -a DIAMOND (Blast default)
rule rgi_run:
    input:
        database = rules.rgi_getdb.output[0],
        input_sequence = filtered_fasta
    output:
        DIR = directory(join(OUTDIR_amr, "rgi/results"))
    params:
        output_file = join(OUTDIR_amr, "rgi/results/rgi"),
        extra = ["--local", "--clean", "--low_quality", 
            "--include_nudge", "--split_prodigal_jobs"]
    log: join(OUTDIR_amr, "rgi/results/rgi.log")
    benchmark: join(OUTDIR_amr, "rgi/results/rgi.bench")
    threads: workflow.cores
    wrapper:
        "file:///local/plsdb/master/wrappers/rgi/run"

# Comparison
#################################################

# The html output seems not to work
rule hamronize:
    input:
        amrfinderplus = rules.AMRFinderPlus_join.output.tsv,
        rgi = join(rules.rgi_run.output.DIR, "rgi.txt"),
    output:
        DIR = directory(join(OUTDIR_amr, "summary")),
        tsv = join(OUTDIR_amr, "combined_report.tsv"),
        html = join(OUTDIR_amr, "combined_report.html")
    params: 
        metadata = config['amr_tools']
    log:
        tsv = join(OUTDIR_amr, "combined_report.log"),
    wrapper:
        "file:///local/plsdb/master/wrappers/hamronization"

rule hamronize_dedup:
    input: rules.hamronize.output.tsv
    # input: "/local/plsdb/test_amr.tsv"
    output: join(OUTDIR_amr, "combined_report_dedup.tsv")
    params:
        range_nt = 10
    log: join(OUTDIR_amr, "combined_report_dedup.log")
    benchmark: join(OUTDIR_amr, "combined_report_dedup.bench")
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/hamronize_derep.py"
    