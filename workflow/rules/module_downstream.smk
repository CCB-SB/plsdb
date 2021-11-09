
# Krona plot
##################################################
rule krona_xml:
    input:
        MASTER_TAB
    output:
        MASTER_KRONA_XML
    message:
        "Create XML file for Krona plot from {input}"
    shell:
        "python scripts/create_krona_xml.py -t {input} -o {output}"

rule krona_html:
    input:
        MASTER_KRONA_XML
    output:
        MASTER_KRONA_HTML
    message:
        "Create Krona plot from {input}"
    params:
        bin=BIN_KRONA_XML
    shell:
        "{params.bin} {input} -o {output}"

# Summary
##################################################
rule summary:
    input:
        tab=MASTER_TAB,
        script='scripts/summary.R'
    output:
        pdf=MASTER_SUM_PLOTS,
        txt=MASTER_SUM_STATS
    message:
        "Summary of collected data: plots and statistics"
    params:
        width=10,
        height=6
    shell:
        "Rscript {input.script} --tab {input.tab} --pdf {output.pdf} --width {params.width} --height {params.height} | tee {output.txt}"


# Compare to older version
##################################################
rule compare:
    input:
        new=MASTER_TAB,
        old=config['old_tab'],
        new_nonfiltered=PLASMIDS_FULL1,
        script='scripts/compare_master_tabs.R'
    output:
        MASTER_CHANGES
    shell:
        "Rscript {input.script} -n {input.new} -o {input.old} -f {input.new_nonfiltered} -t {output} -l {output}.log"

# Similar records
##################################################
rule sim_records:
    input:
        MASTER_MASH_DISTS
    output:
        MASTER_SIM
    message:
        'Lits of similar records pairs from {input}'
    run:
        pairs = set()
        with open(input[0], 'r') as ifile, open(output[0], 'w') as ofile:
            for line in tqdm(ifile):
                sID, qID, dist, pv, sh = line.rstrip('\n').split('\t')
                pair = sorted([sID, qID])
                pair = (pair[0], pair[1])
                # same ID -> skip
                if sID == qID:
                    continue
                # group ID already set -> skip
                if pair in pairs:
                    continue
                else:
                    pairs.add(pair)
                    ofile.write('{}\t{}\n'.format(sID, qID))