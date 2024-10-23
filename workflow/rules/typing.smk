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
        out_file = join(OUTDIR, "typing/mob/typer/report.tsv"),
        biomarker_report_file = join(OUTDIR, "typing/mob/typer/biomarker_report.tsv"),
        mge_report_file = join(OUTDIR, "typing/mob/typer/mge_report.tsv"),
        analysis_dir = directory(join(OUTDIR, "typing/mob/typer"))
    params:
        extra = ['--multi', '--keep_tmp']
    threads: workflow.cores
    log: join(OUTDIR, "typing/mob/typer/results.log")
    benchmark: join(OUTDIR, "typing/mob/typer/results.bench")
    wrapper:
        "file:///local/plsdb/master/wrappers/mob_suite/typer"

rule typing_table:
    input:
        mob = rules.mob_typer.output.out_file,
        pmlst = rules.pmlst_join.output[0],
    output:
        join(OUTDIR, "final", "typing.csv")
    run:
        import pandas as pd
        import numpy as np

        mob = pd.read_table(input.mob, low_memory=False)
        mob['NUCCORE_ACC'] = mob.sample_id.str.split(pat=' ', expand=True).loc[:, [0]]
        mob.drop(columns=["sample_id"], inplace=True)
        
        pmlst = pd.read_table(input.pmlst)
        pmlst['PMLST_sequence_type'] = np.where(np.isnan(pmlst['PMLST_sequence_type']), "", "ST"+pmlst['PMLST_sequence_type'].astype(str))
        df = pd.merge(mob, pmlst, how='left', on="NUCCORE_ACC", validate="1:1")
        df = df.replace({"-":""})
        df = df.fillna(np.nan).replace([np.nan], [None])

        df.to_csv(output[0], index=False)

rule typing_markers_table:
    input:
        rules.mob_typer.output.biomarker_report_file
    output:
        join(OUTDIR, "final", "typing_markers.csv")
    run:
        import pandas as pd
        
        df = pd.read_table(input[0], low_memory=False)
        df['NUCCORE_ACC'] = [i.split(' ')[0] for i in df['sseqid']]
        df['MOB_suite_ID'] = [i.split('|')[0] for i in df['qseqid']]
        df['element'] = [i.split('|')[1] for i in df['qseqid']]

        df.to_csv(output[0], index=False)

# Merge annotations
#################################################

rule features_gbk:
    input:
        fasta = filtered_fasta,
        nucc = filtered_pls,
        bgc = rules.antismash_join.output.tsv,
        genbank = rules.genbank_join.output[0],
        amr = rules.hamronize_dedup.output[0],
        typing = rules.typing_markers_table.output[0]
    output:
        amr_tab = join(OUTDIR, "final", "amr.tsv"),
        gc_tab = join(OUTDIR, "filtering/metadata/nucc_gc.csv"),
        proteins_tab = join(OUTDIR, "final/proteins.csv"),
        proteins = join(OUTDIR, "final/proteins.fasta"),
        DIR = directory(join(OUTDIR, "final/features/gbk/"))
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/features.py"

#git clone https://github.com/stothard-group/cgview-builder.git -b master {output.DIR_repo}
rule features_json:
    input: 
        DIR = rules.features_gbk.output.DIR,
        config = "../src/cgview_config.yaml"
    output:
        DIR = directory(join(OUTDIR, "final/features/json")),
    conda: "../envs/cgview_build.yaml"
    log: join(OUTDIR, "features/features/json.log")
    threads: 120
    params:
        DIR_repo = directory("scripts/cgview-builder/")
    script:
        "../scripts/cgview_json.py"