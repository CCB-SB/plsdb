rule pmlst_split:
    input:
        fastx=filtered_fasta
    output:
        DIR=temp(directory(join(OUTDIR, "pubmlst/pmlst/split/")))
    log:
        join(OUTDIR, "pubmlst/pmlst/split/multi_fasta_split.log"),
    benchmark:
        join(OUTDIR, "pubmlst/pmlst/split/multi_fasta_split.bench"),
    threads: 1
    shell:
        """
        awk '/^>/ {{ file="{output.DIR}/" substr($1,2) ".fasta" }} {{ print > file }}' {input.fastx}
        """
# pMLST
#################################################
rule pmlst_db:
    params:
        db = config['pubmlst']['pmlst']['db'],
        url = config['pubmlst']['pmlst']['url'],
        url_schemes = config['pubmlst']['pmlst']['url_schemes'],
    output: 
        DIR = directory(join(OUTDIR, "typing/pmlst/db/"))
    threads: 4 # maximum for this API
    log: join(OUTDIR, "pmlst", "db/api.log")
    benchmark: join(OUTDIR, "pmlst", "db/api.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/pubmlst/pmlst/db"

# Adapted from https://github.com/tseemann/mlst/blob/master/scripts/mlst-make_blast_db
rule pmlst_blastdb:
    input:
        rules.pmlst_db.output.DIR
    output:
        directory(join(OUTDIR, "typing/pmlst/blast"))
    conda: "../envs/mlst.yaml"
    shell:
        """
        MLSTDIR={input[0]}
        BLASTDIR={output[0]}
        BLASTFILE="$BLASTDIR/mlst.fa"

        mkdir -p "$BLASTDIR"
        rm -f "$BLASTFILE"

        for N in $(find $MLSTDIR -mindepth 1 -maxdepth 1 -type d); do
        SCHEME=$(basename $N)
            echo "Adding: $SCHEME"
            cat "$MLSTDIR"/$SCHEME/*.tfa \
                | grep -v 'not a locus'  \
                | sed -e "s/^>/>$SCHEME./" \
                >> "$BLASTFILE"
        done

        makeblastdb -hash_index -in "$BLASTFILE" -dbtype nucl -title "PubMLST" -parse_seqids

        echo "Created BLAST database for $BLASTFILE"
        """

rule pmlst_run:
    input:
        blastdb = rules.pmlst_blastdb.output[0],
        pmlst_data = rules.pmlst_db.output.DIR,
        fasta_file = join(rules.bgc_split.output.DIR, "{sample}.fasta"),
    output:
        join(OUTDIR, "typing/pmlst/run/{sample}/mlst_res.tsv")
    params:
        threads = 1,
        minid = 85,
        mincov = 66,
        minscore = 50,
        blastdb = lambda wildcards, input: join(input.blastdb, "mlst.fa"),
        datadir = lambda wildcards, input: input.pmlst_data
    threads: 1
    log: join(OUTDIR, "typing/pmlst/run/{sample}/results.log")
    conda: "../envs/mlst.yaml"
    shell:
        """
        mlst --minid {params.minid} --mincov {params.mincov} \
            --minscore {params.minscore} --blastdb {params.blastdb} \
            --datadir {params.datadir} \
            {input.fasta_file} > {output[0]} 2> {log}
        """ 

rule pmlst_join:
    input:
        expand(rules.pmlst_run.output[0], 
            sample=get_nuccore_accs(rules.bgc_split.output.DIR))
    output:
        join(OUTDIR, "typing/pmlst/summary/results.tsv")
    log: 
        join(OUTDIR, "typing/pmlst/summary/results.log")
    run:
        import pandas as pd
        import pdb
        l = []
        
        for f in input:
            d = pd.read_table(str(f), header=None, na_values=['-'])
            d = d.rename({0: "filename", 1: "PMLST_scheme", 2:"PMLST_sequence_type"}, axis=1)
            if d['PMLST_scheme'].dropna().empty:
                continue
            d.insert(loc=0,
                column = "NUCCORE_ACC",
                value=basename(dirname(str(f))) )
            d['PMLST_alleles'] = d.iloc[:, 4:].dropna(axis=1, how='all').agg(','.join,axis="columns")
            d = d[["NUCCORE_ACC","PMLST_scheme", "PMLST_sequence_type", "PMLST_alleles"]]
            l.append(d)

        df = pd.concat(l, ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False)

# MOB
#################################################
rule mob_db:
    output:
        database_directory = directory(join(OUTDIR, "typing/mob/db"))
    log:
        join(OUTDIR, "typing/mob/db/mob_db.log")
    benchmark:
        join(OUTDIR, "typing/mob/db/mob_db.bench")
    threads: 1
    wrapper:
        "file:///local/plsdb/master/wrappers/mob_suite/db"

rule mob_typer:
    input:
        database_directory = rules.mob_db.output.database_directory,
        infile = filtered_fasta
    output:
        out_file = join(OUTDIR, "final/mob.tsv")
    params:
        extra = ['--multi']
    threads: workflow.cores
    log: join(OUTDIR, "typing/mob/typer/results.log")
    benchmark: join(OUTDIR, "typing/mob/typer/results.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/mob_suite/typer"

# CONJSCAN
#################################################
# rule conjscan_model:
#     output:
#         model = directory(join(OUTDIR, "conjscan/model/"))
#     conda: "../envs/macsyfinder.yaml"
#     shell:
#         """
#         macsyfinder install CONJSCAN --models-dir {output.model}
#         """
    
# rule conjscan_run:
#     input:
#         fasta = filtered_fasta,
#         model = rules.output.model
#     output:
#     conda: "../envs/macsyfinder.yaml"
#     shell:
#         """
#         macsyfinder --db-type ordered_replicon \
#                 --sequence-db {input.fasta} \
#                 --models CONJSCAN/Plasmids \
#                 --models-dir {input.model} \
#                 system --option
#         """