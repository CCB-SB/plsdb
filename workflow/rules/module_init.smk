
## abricate: data
rule update_abricate:
    output:
        ABRICATE_UPDATE
    params:
        dbs=config['abricate']['dbs']
    message:
        "Update ABRicate databases: {params.dbs}"
    run:
        logger = setup_logger(logging.INFO)

        for db in params.dbs:
            # run ABRicate to update the databases
            cmd = "abricate-get_db --db {db} --force".format(db=db)
            logger.info('Update with ABRicate DB {}: {}'.format(db, cmd))
            cmd, cmd_s, cmd_o = run_cmd(cmd)
            assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd, cmd_s, cmd_o)
            with open(output[0], 'a') as ofile:
            	ofile.write(cmd_o)

            # reformat FASTA: remove non-UTF-8 characters
            seq0 = "{env_path}/db/{db}/sequence"
            seq1 = seq0 + ".tmp"
            cmd = "mv {input} {input}.tmp && iconv -c -t UTF-8 < {input}.tmp > {input}".format(
                input=os.path.join(ENV_PATH, 'db', db, 'sequences')
            )
            logger.info('Remove non-UTF-8 symbols: {}'.format(cmd))
            cmd, cmd_s, cmd_o = run_cmd(cmd)
            assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd, cmd_s, cmd_o)

        with open(output[0], 'a') as ofile:
            ofile.write('\nUPDATE DONE.\n')

# mlst
# pmlst: data
rule update_pmlst:
    output:
        PMLST_UPDATE
    shell:
        "touch {output}"

rule clean_pmlst:
    input:
        PMLST_UPDATE
    output:
        PMLST_CLEAN
    message:
        "Removing all schemes from mlst db"
    params:
        dir=ENV_PATH
    shell:
        """
        find {params.dir}/db/pubmlst/ -maxdepth 1 -mindepth 1 -type d -exec rm -r {{}} \; && \
        find {params.dir}/db/blast/ -type f -name 'mlst.*' -exec rm {{}} \; && \
        touch {output}
        """