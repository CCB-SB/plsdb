rule bio_table:
    input:
        bio = rules.filter_metadata.output.bio,
        dis_eco = rules.disease_infer_check.params.mapping_checked,
        loc = rules.locations_infer.params.loc_infer,
    output:
        join(OUTDIR, "final", "biosample.csv")
    conda: "../envs/py_env.yaml"
    script:
        "../scripts/biosample_join.py"

rule disease_table:
    input:
        bio = rules.bio_table.output[0],
        do = rules.disease.output.do_terms,
        symp = rules.disease.output.symp_terms,
    output:
        join(OUTDIR, "final", "disease_terms.csv")
    run:
        import pandas as pd

        df = pd.read_csv(input.bio)
        IDs = set(j for i in set(df['DISEASE_ontid']) for j in i.split(',') if i)
        print(f"# of DOID/SYMP: {len(IDs)}")

        # Merge DO and SYMP terms
        df = pd.read_csv(input.symp)
        do = pd.read_csv(input.do)
        if 'DOID:225' in do['id']:
            print("Is present")
        extra_cols = ["derives_from", "has_material_basis_in", "has_symptom"]
        df[extra_cols] = ""
        df = pd.concat([df, pd.read_csv(input.do)])

        # Filter by terms present in dataset
        [i if i in list(df['id']) else print(i) for i in IDs]
        df = df[df['id'].isin (IDs)]
        print(f"Found: {len(df.index)}.Not found printed above")

        df.to_csv(output[0], index=False)


rule taxid_table:
    input:
        tax = rules.filter_metadata.output.tax,
        eco_tax = rules.ecosystem_taxid.output.tax,
    output:
        join(OUTDIR, "final", "taxonomy.csv")
    run:
        import pandas as pd

        df = pd.concat([
            pd.read_csv(input.tax),
            pd.read_csv(input.eco_tax)],
            ignore_index=True)
        
        for cname in df.columns:
            if re.fullmatch(r'TAXONOMY_.*_id', cname, re.IGNORECASE):
                df.loc[:, cname] = df[cname].fillna(-1).astype(int)

        df.drop_duplicates(inplace=True)
        if sum(df.duplicated('TAXONOMY_UID')) > 0:
            print(df[df.duplicated('TAXONOMY_UID', keep=False)])
            raise ValueError
        df.to_csv(output[0], index=False)

rule nucc_table:
    input:
        nucc = filtered_pls,
        gc = rules.features_gbk.output.gc_tab,
        viral = rules.viral_curation.output.pls,
        pmlst = join(OUTDIR, "typing/pmlst/summary/results.tsv") # rules.pmlst_join.output[0],
    output:
        join(OUTDIR, "final", "nuccore.csv")
    run:
        import pandas as pd

        nucc = pd.read_csv(input.nucc)
        nucc.drop(columns=['NUCCORE_Topology'], inplace=True)
        gc = pd.read_csv(input.gc).drop_duplicates()
        viral = pd.read_csv(input.viral)
        pmlst = pd.read_table(input.pmlst)

        print(nucc.columns)
        df1 = pd.merge(nucc, gc, how='left', on="NUCCORE_ACC", validate="1:1")
        print(df1.columns)
        df2 = pd.merge(df1, viral, how='left', on="NUCCORE_ACC", validate="1:1")
        print(df2.columns)
        final = pd.merge(df2, pmlst, how='left', on="NUCCORE_ACC", validate="1:1")
        print(final.columns)

        final.to_csv(output[0], index=False)

rule final_fasta:
    input: filtered_fasta
    output: join(OUTDIR, "final/sequences.fasta")
    run:
        shell("ln -sr {input[0]} {output[0]}")

rule ecopaths_table:
    input: config["ecosystem"]['paths']
    output: join(OUTDIR, "final/ecopaths.csv")
    run:
        shell("ln -sr {input[0]} {output[0]}")

rule createmash:
    input: filtered_fasta
    output:
        DIR = directory(join(OUTDIR, "final/mash")),
        sketch = join(OUTDIR, "final/mash/sketch.msh"),
        distS = join(OUTDIR, "final/mash/mash_distS.tsv"),
        sim = join(OUTDIR, "final/mash/mash_sim.tsv"),
        dist = join(OUTDIR, "final/mash/mash_dist.tsv"),
    conda: "../envs/mash.yaml"
    threads: workflow.cores
    script:
        "../scripts/CreateMashDB.py"


# Summary
##################################################
# rule dstream_summary:
#     input:
#         rules.process_infotable.output
#     output:
#         pdf = join(OUTDIR_dstream, "summary.pdf"),
#         txt = join(OUTDIR_dstream, "summary.txt")
#     params:
#         width=10,
#         height=6
#     conda:
#         "../envs/requirements.yml"
#     shell:
#         """
#         Rscript scripts/dstream_summary.R --tab {input} \
#         --pdf {output.pdf} --width {params.width} \
#         --height {params.height} | tee {output.txt}
#         """


# # Compare to older version
# ##################################################
# rule dstream_compare:
#     input:
#         new = rules.process_infotable.output,
#         old = config['previous_table'],
#         new_nonfiltered = rules.retrieve_plasmid_taxid.output[1]
#     conda:
#         "../envs/requirements.yml"
#     output:
#         txt = join(OUTDIR_dstream, "changes.tsv"),
#         log = join(OUTDIR_dstream, "changes.tsv.log")
#     shell:
#         """
#         Rscript scripts/dstream_compare_tabs.R  -n {input.new} \
#         -o {input.old} -f {input.new_nonfiltered} \
#         -t {output.txt} -l {output.log}
#         """

# # Server data
# ##################################################
# rule dstream_server_data:
#     input:
#         abr = rules.process_join_abricate.output.tsv,
#         changes = rules.dstream_compare.output.txt,
#         changes_log = rules.dstream_compare.output.log,
#         html = rules.dstream_krona_html.output,
#         msh = rules.process_mash_sketch.output,
#         sim = rules.dstream_sim_records.output,
#         infotab = rules.process_infotable.output,
#         fasta = rules.filter_artifacts.output.fasta
#     output:
#         dir = directory("../src/server_data"),
#         fasta = "../src/server_data/plsdb.fna"
#     conda: 
#         "../envs/py_env.yml"
#     shell:
#         """
#         mkdir -p {output.dir} && \
#             cp {input.abr} {output.dir}/plsdb.abr &&\
#             cp {input.changes} {output.dir}/plsdb_changes.tsv &&\
#             cp {input.html} {output.dir}/plsdb.html &&\
#             cp {input.msh} {output.dir}/plsdb.msh &&\
#             cp {input.sim} {output.dir}/plsdb.sim &&\
#             cp {input.infotab} {output.dir}/plsdb.tsv && \
#             cp {input.fasta} {output.dir}/plsdb.fna && \
#             bzip2 -zk {input.fasta} --stdout > {output.dir}/plsdb.fna.bz2
#         """

# # BLAST DBs
# ##################################################
# use rule process_make_rmlst_blastdb as dstream_blastndb with:
#     input:
#         fasta=rules.dstream_server_data.output.fasta,
#     output:
#         dbs = expand(["{file}.{ext}"],
#             file = join("../src/server_data/", f"plsdb.fna"),
#             ext=['nin', 'nhr', 'nsq', 'ndb', 'njs', 'not', 'ntf', 'nto'])
#     params:
#         title = 'plsdb'