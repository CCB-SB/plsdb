OUTDIR_bgc = join(OUTDIR, "bgc")

rule bgc_split:
    input:
        fastx=filtered_fasta
    output:
        DIR=directory(join(OUTDIR_bgc, "split/"))
    log:
        join(OUTDIR_bgc, "split/multi_fasta_split.log"),
    benchmark:
        join(OUTDIR_bgc, "split/multi_fasta_split.bench"),
    threads: 1
    shell:
        """
        awk '/^>/ {{ file="{output.DIR}/" substr($1,2) ".fasta" }} {{ print > file }}' {input.fastx}
        """
def get_nuccore_accs(DIR):
    from glob import glob
    nuccores = [splitext(basename(i))[0] for i in glob(join(DIR, "*.fasta"))]
    return nuccores


# ANTISMASH
#################################################
rule antismash_docker_script:
    output: "scripts/bgc/antismash_run.sh"
    shell:
        """
        curl -q https://dl.secondarymetabolites.org/releases/latest/docker-run_antismash-full > {output} && \
        chmod a+x {output}
        """

rule antismash_run:
    input:
        antismash_script = rules.antismash_docker_script.output,
        fna = join(rules.bgc_split.output.DIR, "{part}.fasta")
    output:
        gbk = join(OUTDIR, 'final/antismash/{part}/{part}.gbk'),
        json = join(OUTDIR, 'final/antismash/{part}/{part}.json'),
        DIR = directory(join(OUTDIR, 'final/antismash/{part}/'))
    params:
        outdir = lambda wildcards, output: dirname(output.DIR),
        basename = lambda wildcards: wildcards.part
    log: join(OUTDIR_bgc, 'antismash/{part}.log')
    threads: 1
    wrapper:
        "file:///local/plsdb/master/wrappers/antismash/run"  

rule antismash_summary:
    input:
        lambda wildcards: 
            expand(rules.antismash_run.output.DIR, part=wildcards.part)
    output:
        tsv = join(OUTDIR, "final/antismash/summary", "{part}.tsv")
    threads: 1
    wrapper:
        "file:///local/plsdb/master/wrappers/antismash/summary"

ruleorder: antismash_run > antismash_summary
rule antismash_join:
    input:
        expand(rules.antismash_summary.output.tsv,
            part=get_nuccore_accs(rules.bgc_split.output.DIR))
    output:
        tsv = join(OUTDIR, "final/antismash/summary", "plsdb.tsv")
    threads: 1
    run:
        import pandas as pd
        df = pd.concat([pd.read_table(str(f)) for f in input ], ignore_index=True)
        df.to_csv(output.tsv, sep="\t", index=False)


# # BIGSCAPE
# #################################################

# rule bigscape_download_pfam_db:
#     output:
#         dir = directory(join(OUTDIR, "bigscape/pfam/"))
#     params:
#         link = config['pfam']['release']
#     conda: "../envs/bgc/bigscape.yml"
#     shell:
#         """
#         mkdir {output.dir}
#         wget -P {output.dir} {params.link}
#         gunzip -c {output.dir}/Pfam-A.hmm.gz > {output.dir}/Pfam-A.hmm
#         hmmpress {output.dir}/Pfam-A.hmm
#         """

# # NOTE: current version of conda has a bug: https://github.com/medema-group/BiG-SCAPE/issues/90
# rule bigscape_download:
#     output: 
#         dir_source = directory(join(OUTDIR, "bigscape", "source", 
#             f"BiG-SCAPE-{config['bigscape']['version']}"))
#     params:
#         zip_path = join(OUTDIR, "bigscape", "source"),
#         link = config['bigscape']['release']
#     shell:
#         """
#         wget {params.link} -O {params.zip_path}.zip
#         unzip {params.zip_path}.zip -d {params.zip_path}
#         rm {params.zip_path}.zip
#         """

# rule bigscape:
#     input:
#         script = rules.bigscape_download.output.dir_source,
#         input_dir = rules.antismash.output.dir,
#         pfam_dir = rules.bigscape_download_pfam_db.output.dir
#     output:
#         dir = join(OUTDIR, "bigscape/results")
#     params:
#         cutoffs = config['bigscape']['cutoffs'],
#         clan_cutoff = config['bigscape']['clan_cutoff']
#     benchmark:
#          join(OUTDIR, "bigscape", "benchmark.txt")
#     log: join(OUTDIR, "bigscape", "log.txt")
#     threads:config['bigscape']['threads']
#     conda: "../envs/bgc/bigscape.yml"
#     shell:
#         """
#         python {input.script}/bigscape.py --verbose --mibig \
#             --include_gbk_str region HMBGC \
#             --cutoffs {params.cutoffs} \
#             --clan_cutoff {params.clan_cutoff} \
#             --pfam_dir {input.pfam_dir} \
#             --cores {threads} \
#             --inputdir {input.input_dir} \
#             --outputdir {output.dir} \
#             &> "{log}"
#         """

# rule bigscape_summary:
#     input: 
#         DIR = rules.bigscape.output.dir
#     output: 
#         tsv = join(OUTDIR, "summary", "bigscape_summary.tsv")
#     params:
#         cutoff_gcf=config['bigscape']['clan_cutoff'][0],
#         cutoff_gcc=config['bigscape']['clan_cutoff'][1]
#     conda: "../envs/bgc/antismash.yml"
#     script:
#         "../scripts/bigscape_summary.py"
 
