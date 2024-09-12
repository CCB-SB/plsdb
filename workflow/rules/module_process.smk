
# NUCCORE GC content
##################################################
rule process_calculate_GC:
    input:
        pls = rules.filter_artifacts.output.md,
        fasta = rules.filter_artifacts.output.fasta
    output:
        join(OUTDIR_process, "calculate_gc.csv")
    log:
        join(OUTDIR_logs, "process_calculate_gc.log")
    conda:
        "../envs/py_env.yml"
    script:
        "../scripts/process_calculate_gc.py"

# UMAP
#################################################

# # Embedding using UMAP on Mash distances
# rule process_umap:
#     input:
#         rules.process_mash_dist.output
#     output:
#         join(OUTDIR_process, "umap_mash_dist.umap")
#     params:
#         neighbors=config['umap']['neighbors'],
#         components=config['umap']['components'],
#         min_dist=config['umap']['min_dist']
#     log:
#         join(OUTDIR_logs, "process_umap.log")
#     conda:
#         "../envs/requirements.yml"
#     shell:
#         """
#         python scripts/process_umap.py \
#             --input {input} --ofile {output} \
#             --log {log} --neighbors {params.neighbors} \
#             --components {params.components} \
#             --min_dist {params.min_dist}
#         """
        

# # Info table
# ##################################################
# rule process_infotable:
#     input:
#         pls = rules.process_parse_locations.output.pls,
#         location_checking = rules.process_manually_inspect_locations.output,
#         emb = rules.process_umap.output,
#         abr = rules.process_join_abricate.output.csv,
#         pmlst = rules.process_pmlst.output
#     output:
#         join(OUTDIR_process, f"infotable_{config['version']}.tsv")
#     log:
#         join(OUTDIR_logs, "process_infotable.log")
#     conda:
#         "../envs/py_env.yml"
#     shell:
#         """
#         python scripts/process_infotable.py --pls {input.pls} \
#             --emb {input.emb} --abr {input.abr} \
#             --pmlst {input.pmlst} --ofile {output} \
#             --log {log}
#         """
