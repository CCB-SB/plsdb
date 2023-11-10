
# Krona plot
##################################################
# "Create XML file for Krona plot from {input}"
rule dstream_krona_xml:
    input:
        rules.process_infotable.output
    output:
        join(OUTDIR_dstream, "krona.xml")
    log:
        join(OUTDIR_logs, "dstream_krona_xml.log")
    conda:
        "../envs/requirements.yml"
    shell:
        "python3 scripts/dstream_krona_xml.py -t {input} -o {output} --log {log}"

# "Create Krona plot from {input}"
rule dstream_krona_html:
    input:
        rules.dstream_krona_xml.output
    output:
        join(OUTDIR_dstream, "krona.html")
    log:
        join(OUTDIR_logs, "dstream_krona_html.log")
    conda:
        "../envs/requirements.yml"
    shell:
        "ktImportXML {input} -o {output} 2> {log}"


# Summary
##################################################
rule dstream_summary:
    input:
        rules.process_infotable.output
    output:
        pdf = join(OUTDIR_dstream, "summary.pdf"),
        txt = join(OUTDIR_dstream, "summary.txt")
    params:
        width=10,
        height=6
    conda:
        "../envs/requirements.yml"
    shell:
        """
        Rscript scripts/dstream_summary.R --tab {input} \
        --pdf {output.pdf} --width {params.width} \
        --height {params.height} | tee {output.txt}
        """


# Compare to older version
##################################################
rule dstream_compare:
    input:
        new = rules.process_infotable.output,
        old = config['previous_table'],
        new_nonfiltered = rules.retrieve_plasmid_taxid.output[1]
    conda:
        "../envs/requirements.yml"
    output:
        txt = join(OUTDIR_dstream, "changes.tsv"),
        log = join(OUTDIR_dstream, "changes.tsv.log")
    shell:
        """
        Rscript scripts/dstream_compare_tabs.R  -n {input.new} \
        -o {input.old} -f {input.new_nonfiltered} \
        -t {output} -l {output}.log
        """


# Similar records
##################################################
rule dstream_sim_records:
    input:
        rules.process_mash_dist_sim.output
    output:
        join(OUTDIR_process, "mash", "mash_dist.sim")
    conda:
        "../envs/requirements.yml"
    shell:
        "python scripts/dstream_sim_records.py --input {input} --ofile {output}"

# BLAST DBs
##################################################
use rule process_make_rmlst_blastdb as dstream_blastndb with:
    input:
        fasta=rules.filter_artifacts.output.fasta,
    output:
        dbs = expand(["{file}.{ext}"],
            file = join(OUTDIR_filtering, f"artifact_filtering.fna"),
            ext=config['blast']['ext'])
    params:
        title = f"plsdb_{config['version']}"


# Server data
##################################################
rule dstream_server_data:
    input:
        abr = rules.process_join_abricate.output.tsv,
        changes = rules.dstream_compare.output.txt,
        changes_log = rules.dstream_compare.output.log,
        nin = rules.dstream_blastndb.output[0],
        nhr = rules.dstream_blastndb.output[1],
        nsq = rules.dstream_blastndb.output[2],
        html = rules.dstream_krona_html.output,
        msh = rules.process_mash_sketch.output,
        sim = rules.dstream_sim_records.output,
        infotab = rules.process_infotable.output,
        fasta = rules.filter_artifacts.output.fasta
    output:
        dir = directory("../src/server_data")
    conda: 
        "../envs/py_env.yml"
    shell:
        """
        mkdir -p {output.dir} && \
            cp {input.abr} {output.dir}/plsdb.abr &&\
            cp {input.changes} {output.dir}/plsdb_changes.tsv &&\
            cp {input.nin} {output.dir}/plsdb.fna.nin &&\
            cp {input.nhr} {output.dir}/plsdb.fna.nhr &&\
            cp {input.nsq} {output.dir}/plsdb.fna.nsq &&\
            cp {input.html} {output.dir}/plsdb.html &&\
            cp {input.msh} {output.dir}/plsdb.msh &&\
            cp {input.sim} {output.dir}/plsdb.sim &&\
            cp {input.infotab} {output.dir}/plsdb.tsv && \
            bzip2 -zk {input.fasta} --stdout > {output.dir}/plsdb.bz2
        """