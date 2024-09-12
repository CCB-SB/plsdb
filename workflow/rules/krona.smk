# Krona plot
##################################################
# NOTE: for some reason conda env is not installing krona. For now, manual installation in env
rule krona_ncbidb:
    output: 
        DIR = directory(join(OUTDIR, "krona/ncbidb/"))
    conda: "../envs/krona.yml"
    log: join(OUTDIR, "krona/ncbidb/krona_ncbidb.log")
    shell:
        """
        ktUpdateTaxonomy.sh {output.DIR} 2> {log}
        """    
# NOTE: need to pass to tsv
rule krona_taxonomy:
    input:
        db = rules.krona_ncbidb.output.DIR,
        tax_tsv = "../../results/filtering/metadata/tax_filt.tsv" #rules.filter_metadata.output.tax
    output:
        html = join(OUTDIR, "final/krona.html")
    params:
        taxonomy_uid_col = 1, # TAXONOMY_UID
        query_id_col = 20, # TAXONOMY_taxon_lineage
    log:
        join(OUTDIR, "krona/taxonomy/taxonomy_krona.log")
    conda: "../envs/krona.yml"
    shell:
        """
        sed '1d' {input.tax_tsv} > tmp.taxonomy.tsv # Remove header
        ktImportTaxonomy tmp.taxonomy.tsv -tax {input.db} \
            -t {params.taxonomy_uid_col} -q {params.query_id_col} \
            -o {output.html} 2> {log}
        rm tmp.taxonomy.tsv
        """
