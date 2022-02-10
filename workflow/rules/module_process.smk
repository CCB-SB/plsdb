
rule create_pmlst_db:
    input:
        PMLST_ALLELES
    output:
        PMLST_DB
    message:
        "Creating pMLST db"
    params:
        dir=ENV_PATH
    shell:
        "{params.dir}/scripts/mlst-make_blast_db"

# rmlst: data
# NOTE: FASTA file must exist
rule rmlst_blastdb:
    input:
        fasta=RMLST_FASTA,
        bin=BIN_MAKEBLASTDB
    output:
        RMLST_DB
    message:
        "Create BLASTDB from {input}"
    shell:
        "{input.bin} -in {input.fasta} -input_type fasta -dbtype nucl -title \"rmlst\""


rule proc_taxon:
    input:
        PLASMIDS_TAX
    output:
        PLASMIDS_TAXP
    message:
        "Process taxa from {input}"
    run:
        logger = setup_logger(logging.INFO)

        header = config['eutils']['header']['taxon']
        ranks  = config['eutils']['header']['ranks']
        for rank in ranks:
            header += ['taxon_%s_id' % rank, 'taxon_%s_name' % rank]

        df = pandas.DataFrame(columns=header)
        with open(input[0], 'r') as ifile:
            i = 0
            for line in tqdm(ifile):
                qids  = []
                qinfo = None
                line = line.rstrip('\n')
                for field in line.split('\t'):
                    if re.search('\|',field):
                        fields = field.split('|')

                        # taxon ID, name, rank
                        tid   = None
                        tname = fields[-2]
                        trank = fields[-1]
                        if trank == 'no rank':
                            trank = 'NA'

                        # taxon ID
                        if len(qids) == 0: # query ID, save all
                            for f_i, f_n in enumerate(fields):
                                if re.fullmatch(r'\d+', f_n) and f_i < (len(fields) - 2):
                                    qids.append(f_n)
                            if len(qids) > 1:
                                logger.info('Taxon with multiple IDs: {}'.format(field))
                            qinfo = dict.fromkeys(qids)
                            for qid in qids:
                                qinfo[qid] = dict.fromkeys(header, 'NA')
                                qinfo[qid]['taxon_id'] = qid
                                qinfo[qid]['taxon_name'] = tname
                                qinfo[qid]['taxon_rank'] = trank
                        else: # ignore mult. IDs, take first one
                            tid = fields[0]

                        # if relevant rank
                        assert qinfo is not None
                        if 'taxon_%s_name' % trank in header:
                            for qid in qids:
                                if tid is None:
                                    qinfo[qid]['taxon_%s_id' % trank] = qid
                                else:
                                    qinfo[qid]['taxon_%s_id' % trank] = tid
                                qinfo[qid]['taxon_%s_name' % trank] = tname
                                qinfo[qid]['taxon_%s_rank' % trank] = trank
                    else: # lineage
                        assert qinfo is not None
                        for qid in qids:
                            qinfo[qid]['lineage'] = field
                # save info for each query
                for qid in qids:
                    df.loc[qid] = [qinfo[qid][k] for k in header]
            # replace nulls by 'NA'
            df.fillna(value='NA', inplace=True)
            # keep only unqiue rows
            df.drop_duplicates(inplace=True)
            # save
            df.to_csv(output[0], sep='\t', header=True, index=False, index_label=False)

# Add to table
rule collect_meta:
    input:
        pls=PLASMIDS_ALL,
        asm=PLASMIDS_ASM,
        bios=PLASMIDS_BIOS,
        tax=PLASMIDS_TAXP,
        pls_asm=PLASMIDS_ALL_ASM,
        pls_bios=PLASMIDS_ALL_BIOS,
    output:
        PLASMIDS_FULL1
    message:
        "Adding meta info to {input.pls}"
    log:
        os.path.join(config['logs']['odir']['reports'], '%s__full1.log' % timestamp)
    run:
        import datetime
        logger = setup_logger(logging.INFO,log[0])
        logger = setup_logger(logging.INFO)
        pls = pandas.read_csv(input.pls, sep='\t', header=0, dtype=str)
        pls.set_index('UID_NUCCORE', drop=False, inplace=True, verify_integrity=True)
        logger.info('Read in {} plasmid records\n{}'.format(pls.shape[0], pls.head()))

        asm = pandas.read_csv(input.asm, sep='\t', header=0, dtype=str)
        asm.set_index(asm.UID_ASSEMBLY, drop=True, inplace=True, verify_integrity=True)
        logger.info('Read in {} assembly records\n{}'.format(asm.shape[0], asm.head()))

        tax = pandas.read_csv(input.tax, sep='\t', header=0, dtype=str)
        tax.set_index('taxon_id', drop=True, inplace=True, verify_integrity=True)
        logger.info('Read in {} taxonomy records\n{}'.format(tax.shape[0], tax.head()))

        bios = pandas.read_csv(input.bios, sep='\t', header=0, dtype=str)
        bios.set_index('UID_BIOSAMPLE', drop=True, inplace=True, verify_integrity=True)
        logger.info('Read in {} biosample records\n{}'.format(bios.shape[0], bios.head()))

        pls_asm = pandas.read_csv(input.pls_asm, sep='\t', header=0, dtype=str)
        pls_asm.set_index('UID_NUCCORE', drop=True, inplace=True, verify_integrity=True)
        logger.info('Read in {} plasmid/assembly records\n{}'.format(pls_asm.shape[0], pls_asm.head()))

        pls_bios = pandas.read_csv(input.pls_bios, sep='\t', header=0, dtype=str)
        pls_bios.set_index('UID_NUCCORE', drop=True, inplace=True, verify_integrity=True)
        logger.info('Read in {} plasmid/assembly records\n{}'.format(pls_bios.shape[0], pls_bios.head()))
        
        # process plasmid/assembly mapping - only latest or none
        def get_latest_assembly(uids):
            if pandas.isnull(uids):
                return None
            uids = uids.split(';')
            if len(uids) == 1:
                assert uids[0] in asm.index, 'Unknown UID, check assembly.txt and nuccore_assembly.txt and consider removal and rerun for assembly.txt: {}'.format(uids[0])
                return uids[0]
            for uid in uids:
                assert uid in asm.index, 'Unknown UID, check assembly.txt and nuccore_assembly.txt and consider removal and rerun for assembly.txt: {}'.format(uid)
                if asm.loc[uid,'Latest_ASSEMBLY'] == 'True':
                    return uid
            # get dates, parse, sort, get first
            dates = asm.loc[uids,'SeqReleaseDate_ASSEMBLY'].map(lambda x: datetime.datetime.strptime(x, '%Y/%m/%d %H:%M'))
            dates.sort_values(ascending=False, inplace=True)
            logger.warn('No latest assembly among {}: latest w.r.t. date is {}'.format(uids, dates.index[0]))
            return dates.index[0]
        #for uid in pls_asm['UID_ASSEMBLY']:
        #    assert uid in asm['UID_ASSEMBLY'], 'Unknown UID {}'.format(uid)
        pls_asm['UID_ASSEMBLY'] = pls_asm['UID_ASSEMBLY'].map(get_latest_assembly)
        # add assembly ID
        pls = pandas.merge(
            left=pls,
            right=pls_asm,
            how='left',
            left_index=True,
            right_index=True,
            sort=False,
        )

        # add assembly info
        pls = pandas.merge(
            left=pls,
            right=asm,
            how='left',
            left_on='UID_ASSEMBLY',
            right_index=True,
            sort=False,
        )

        # add biosample ID
        pls = pandas.merge(
            left=pls,
            right=pls_bios,
            how='left',
            left_index=True,
            right_index=True,
            sort=False,
        )

        # add biosample info
        pls = pandas.merge(
            left=pls,
            right=bios,
            how='left',
            left_on='UID_BIOSAMPLE',
            right_index=True,
            sort=False,
        )

        # add taxonomy
        pls = pandas.merge(
            left=pls,
            right=tax,
            how='left',
            left_on='TaxonID_NUCCORE',
            right_index=True,
            sort=False,
        )

        # check if this column name is really used
        if 'UID_ASSEMBLY_x' in pls:
            # check if columns are really identical
            assert pls['UID_ASSEMBLY_x'].equals(pls['UID_ASSEMBLY_y']), 'Unknown UID contained'
            # merge the two identical columns (else pipe gets problems later on)
            pls = pls.rename(columns={"UID_ASSEMBLY_x": "UID_ASSEMBLY"})
            pls = pls.drop(columns="UID_ASSEMBLY_y")        

        # save
        pls.to_csv(output[0], sep='\t', header=True, index=False, index_label=False)


#TODO fix logfile
# rMLST
##################################################
# run rMLST
rule rmlst_blastn:
    input:
        fasta="{basename}.fna",
        db=RMLST_FASTA,
        dbs=RMLST_DB
    output:
        "{basename}.rmlst"
    params:
        bin=BIN_BLASTN,
        cores=30
#    log:
#        os.path.join(config['logs']['odir']['reports'], '%s__full1.log' % timestamp)
    run:
        logger = setup_logger(logging.INFO)

        header = config['rmlst']['rmlst_header']
        ident  = config['rmlst']['rmlst_ident']
        cov    = config['rmlst']['rmlst_cov']

        # blastn search
        cmd = config['rmlst']['rmlst_cmd'].format(
            blastn=params.bin,
            input=input.fasta,
            db=input.db,
            output=output[0],
            ident=ident,
            header=' '.join(header),
            cores=params.cores
        )
        logger.info(cmd)
        cmd, cmd_s, cmd_o = run_cmd(cmd)
        assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd, cmd_s, cmd_o)

        # read in results
        df = pandas.read_csv(output[0], sep='\t', header=None, names=header)

        # compute coverage
        df['cov'] =  100 * (df['length'] - df['gaps']) / df['slen']
        # filter by coverage
        df = df.loc[df['cov'] >= cov,:]

        # rMLST locus
        df['slocus'] = df['sseqid'].map(lambda x: x.split('_')[0])

        # save
        df.to_csv(
            path_or_buf=output[0],
            sep='\t',
            header=True,
            index=False,
            index_label=False
        )


######################################################################
# MASTER FILES
######################################################################
rule master_fasta:
    input:
        pls=PLASMIDS_FILT4,
        fna=PLASMIDS_FILT3_FASTA
    output:
        MASTER_FASTA
    log:
        os.path.join(config["logs"]["odir"]["master"], "%s.fna.log" % timestamp)
    message:
        "Create FASTA file for records in {input.pls}"
    shell:
        """
        python scripts/subsample_fasta.py -i {input.pls} -f {input.fna} -o {output} -l {log}
        """
        # from Bio import SeqIO

        # logger = setup_logger(logging.INFO)

        # mkdir(os.path.dirname(output[0]), True)

        # ids = list(pandas.read_csv(input.pls, sep='\t', header=0, dtype=str)['ACC_NUCCORE'])
        # n   = len(ids)
        # ids = set(ids)
        # assert len(ids) == n, 'FASTA IDs in {} are not unique'.format(input.pls)
        # logger.info('Read in {} plasmid records'.format(len(ids)))

        # with open(input.fna, 'r') as ifile, open(output[0], 'w') as ofile:
        #     for record in tqdm(SeqIO.parse(ifile, 'fasta')):
        #         if record.id in ids:
        #             SeqIO.write(record, ofile, 'fasta')
        #             ids.remove(record.id)
        # assert len(ids) == 0, 'Not all IDs found in FASTA {}: {}'.format(input.fna, ';'.join(ids))

# ABRicate annot
##################################################
# NOTE: Using a procedure analoguous to that of PlasmidFinder (https://bitbucket.org/genomicepidemiology/plasmidfinder)
#   1) Call Blaster (from CGE core module): call BLAST and pre-process hits (best hit per subject)
#   2) Filter results by identity and coverage; from overlapping hits keep only the best hit
#   Differences/modifications: see tag "CHANGED"
rule abricate:
    input:
        fna="{basename}.fna",
        abr=ABRICATE_UPDATE,
        blastn=BIN_BLASTN
    output:
        "{basename}.abr.{rdb}"
    message:
        "Search {wildcards.rdb} in fasta: {input.fna}"
    params:
        cov=lambda wildcards: config['abricate']['params'][wildcards.rdb]['cov'],
        ident=lambda wildcards: config['abricate']['params'][wildcards.rdb]['ident']
    run:
        from scripts.utils import call_blaster
        logger = setup_logger(logging.INFO)

        db_path = os.path.join(ENV_PATH, 'db', wildcards.rdb)
        logger.info('Blaster database in {}'.format(db_path))

        # call Blaster to get and pre-process BLAST hits
        df_all = []
        df = []
        
        for record in tqdm(SeqIO.parse(input.fna, 'fasta')):
            # call Blaster: get BLAST hits and pre-process them
            result = call_blaster(record=record, fasta=input.fna, db_path=db_path, blast=input.blastn, cov=params.cov, ident=params.ident)
            
            # no results
            if not result:
                continue

            # convert to table
            result = pandas.DataFrame(result)
            
            # filter by coverage and identity
            result = result.loc[(result['cov'] >= params.cov) & (result['pident'] >= params.ident),:]
            # no results left
            if result.shape[0] == 0:
                continue

            # compute score from identity and coverage
            result = result.assign(ic_score = result['pident'] * result['cov'])

            # sort by query start: small -> large (ascending)
            result.sort_values(by=['qstart'], ascending=[True], inplace=True)
            df_all.append(result)

            # overlapping hits
            keep_hit = [True] * result.shape[0] # True = keep, False = remove
            current_i   = 0
            current_end = result['qend'].values[0]
            current_ics = result['ic_score'].values[0]
            for i in range(0, result.shape[0]-1):
                next_start = result['qstart'].values[i+1]
                next_end   = result['qend'].values[i+1]
                next_ics   = result['ic_score'].values[i+1]

                if next_start <= current_end: # overlapping hits; CHANGED: <=, NOT <
                    # remove hit with lower identity-coverage score
                    if current_ics < next_ics:
                        keep_hit[current_i] = False
                        # update current
                        current_i   = i+1
                        current_end = next_end
                        current_ics = next_ics
                    else:
                        keep_hit[i+1] = False
                else: # no overlap, update current; CHANGED: otherwise if no overlap always comparing to the latest "current"
                    current_i   = i+1
                    current_end = next_end
                    current_ics = next_ics
            result = result.loc[keep_hit,:]

            # rm coverage-identity score
            result = result.drop(['ic_score'], axis=1, inplace=False) # inplace=True produces SettingWithCopyWarning

            # save
            df.append(result)
        
        df_all = pandas.concat(df_all, axis=0, sort=False)
        df = pandas.concat(df, axis=0, sort=False)

        # stitle -> DB + ID
        # title format: <BLAST DB ID> <DB>~~~<name>~~~<ID> <name>
        # example: gnl|BL_ORD_ID|578 argannot~~~(Bla)blaFRI-3~~~KY524440:1-885 (Bla)blaFRI-3
        df_all['sseqdb'] = df_all['stitle'].map(lambda x: x.split(' ')[1].split('~~~')[0])
        df_all['sseqid'] = df_all['stitle'].map(lambda x: ', '.join(x.split(' ')[1].split('~~~')[1:]))
        df['sseqdb']     = df['stitle'].map(lambda x: x.split(' ')[1].split('~~~')[0])
        df['sseqid']     = df['stitle'].map(lambda x: ', '.join(x.split(' ')[1].split('~~~')[1:]))

        # save
        df_all.to_csv(
            path_or_buf=output[0]+'.all',
            sep='\t',
            header=True,
            index=False,
            index_label=False
        )
        df.to_csv(
            path_or_buf=output[0],
            sep='\t',
            header=True,
            index=False,
            index_label=False
        )

rule cat_abricate:
    input:
        ["{basename}.abr.%s" % db for db in config['abricate']['dbs']]
    output:
        "{basename}.abr"
    message:
        "Concat ABRicate hits: {input}"
    params:
        blastn=BIN_BLASTN,
        cores=10
    run:
        dfs = []
        for ifile in input:
            dfs.append(pandas.read_csv(ifile, sep='\t', header=0))
        dfs = pandas.concat(dfs, axis=0, sort=False)
        dfs.to_csv(
            path_or_buf=output[0],
            sep='\t',
            header=True,
            index=False,
            index_label=False
        )

# pMLST
#################################################
rule pmlst:
    input:
        db=PMLST_DB,
        fasta="{path}/{date}.fna",
        pf="{path}/{date}.abr"
    output:
        "{path}/{date}.pmlst"
    message:
        "pmlst on fasta: {input}"
    run:
        import shutil
        from Bio import SeqIO
        from scripts.utils import process_pmlst_hits

        logger = setup_logger(logging.INFO)

        # PlasmidFinder
        pf = pandas.read_csv(input.pf, sep='\t', header=0)
        pf = pf.loc[pf['sseqdb'] == 'plasmidfinder',]
        logger.info('Read in {} PlasmidFinder hits\n{}'.format(pf.shape[0], pf.head()))

        # pMLST
        pmlst_hits = []
        tmp_file = 'pmlst.tmp'
        tmp_fasta = os.path.join(tmp_file + '.fasta')
        logger.info('Run pMLST for each record with available scheme w.r.t. PlasmidFinder hits.')
        with open(input.fasta, 'r') as ifile:
            for record in tqdm(SeqIO.parse(ifile, 'fasta')):
                
                # found replicons
                record_replicons = list(pf.loc[pf['qseqid'] == record.id,'sseqid'])
                if len(record_replicons) == 0:
                    continue
                else:
                    SeqIO.write(record, open(tmp_fasta, 'w'), 'fasta') # create tmp FASTA

                used_schemes = set()
                for record_replicon in record_replicons:
                    # set scheme
                    scheme = None
                    for k,v in config['pmlst']['map'].items():
                        if re.match(k, record_replicon, re.IGNORECASE):
                            scheme = v
                            break
                    if scheme is None or scheme in used_schemes: # no matching scheme/already done
                        continue
                    else:
                        used_schemes.add(scheme)
                    
                    # run pmlst
                    cmd = config['pmlst']['cmd'].format(scheme=scheme) + " {fasta} > {ofile}".format(fasta=tmp_fasta, ofile=tmp_file)
                    cmd, cmd_s, cmd_o = run_cmd(cmd)
                    assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd, cmd_s, cmd_o)
                    
                    # process hits
                    hits = process_pmlst_hits(f=tmp_file, pmlst_db_path=os.path.join(ENV_PATH, 'db/pubmlst'), logger=logger)
                    for i in range(0, len(hits)):
                        hits[i] = {
                            'ID': record.id,
                            'pmlst': hits[i]
                        }
                    pmlst_hits += hits
        
        # add pMLST hits to PlasmidFinder hits
        pmlst_hits = pandas.DataFrame(pmlst_hits)
        # pmlst_hits.set_index(keys='ID', drop=True, inplace=True, verify_integrity=True)
        # pfs = pandas.merge(
        #     left=pfs,
        #     right=pmlst_hits,
        #     how='left',
        #     left_index=True,
        #     right_index=True,
        #     sort=False,
        # )

        # save
        pmlst_hits.to_csv(
            path_or_buf=output[0],
            sep='\t',
            header=True,
            index=False,
            index_label=False
        )

        # clean up
        os.remove(tmp_file)
        os.remove(tmp_fasta)

# MOB_suite
#################################################
rule mob_typer:
    input:
        fasta="{path}/{date}.fna"
    output:
        "{path}/{date}.mob"
    message:
        "mobtyper on fasta: {input}"
    conda:
        "../envs/mobtyper.yml"
    threads:
        64
    shell:
        "mob_typer --num_threads {threads} --multi --infile {input.fasta} --out_file {output}"

# BLASTn DBs
##################################################
# See "general" rules

# Mash
##################################################
# See "general" rules



# Embedding
##################################################
# See "general" rules

# Info table
##################################################
rule infotable:
    input:
        pls=PLASMIDS_DIS_ONT,
        emb=MASTER_MASH_UMAP,
        abr=MASTER_ABRICATE,
        pmlst=MASTER_PMLST,
        mob=MASTER_MOB
    output:
        MASTER_TAB
    message:
        "Create info table for records in {input.pls}"
    run:
        import re
        logger = setup_logger(logging.INFO)

        # plasmids
        pls = pandas.read_csv(input.pls, sep='\t', header=0, dtype=str)
        pls.set_index(keys='ACC_NUCCORE', drop=False, inplace=True, verify_integrity=True)
        logger.info('Read in {} plasmid records\n{}'.format(pls.shape[0], pls.head()))

        # embedding
        emb = pandas.read_csv(input.emb, sep='\t', header=0, index_col='ID', dtype=str)
        logger.info('Read in embedding for {} records\n{}'.format(emb.shape[0], emb.head()))

        # PlasmidFinder
        abr = pandas.read_csv(input.abr, sep='\t', header=0)
        abr = abr.loc[abr['sseqdb'] == 'plasmidfinder',['qseqid', 'sseqid']]
        logger.info('Read in {} PlasmidFinder hits\n{}'.format(abr.shape[0], abr.head()))
        pf = []; pf_index = []
        for aggr_id, aggr_df in abr.groupby(['qseqid']):
            assert len(set(aggr_df['sseqid'])) == len(aggr_df['sseqid'])
            pf_index.append(list(aggr_df['qseqid'])[0])
            pf.append({'plasmidfinder': '|'.join(aggr_df['sseqid'])})
        pf = pandas.DataFrame(pf, index=pf_index)
        logger.info('Aggr PlasmidFinder hits\n{}'.format(pf.head()))

        # pMLST
        pmlst_df = pandas.read_csv(input.pmlst, sep='\t', header=0)
        logger.info('Read in {} pMLST records\n{}'.format(pmlst_df.shape[0], pmlst_df.head()))
        pmlst = []; pmlst_index = []
        for aggr_id, aggr_df in pmlst_df.groupby(['ID']):
            assert len(set(aggr_df['pmlst'])) == len(aggr_df['pmlst'])
            pmlst_index.append(list(aggr_df['ID'])[0])
            pmlst.append({'pmlst': '|'.join(aggr_df['pmlst'])})
        pmlst = pandas.DataFrame(pmlst, index=pmlst_index)
        logger.info('Aggr pMLST hits\n{}'.format(pmlst.head()))

        #mob
        mob = pandas.read_csv(input.mob, sep='\t', header=0)[["sample_id","relaxase_type(s)","mpf_type"]]
        #regex=re.compile('^([A-Z][A-Z]_)?([A-Z][A-Z])?\d+\.\d')

        def extract_id(x):
                regex=re.compile('\d+\.\d$')
                tmp=x.split("_")
                result=tmp[0]
                counter=1
                while not regex.search(result):
                     result=result+"_"+tmp[counter]
                     counter=counter+1	
                return(result)

        mob['sample_id'] = mob['sample_id'].apply(lambda x: extract_id(x) )
        mob[mob=="-"]=""
        mob=mob.set_index("sample_id")
        logger.info('Read in {} Mobtyper records\n{}'.format(mob.shape[0], mob.head()))

        # process source strings
        source_map = {
            'insd': 'INSDC',
            'refseq': 'RefSeq'
        }
        pls['Source_NUCCORE'] = pls['Source_NUCCORE'].map(lambda x: source_map[x])

        # add data
        pls = pandas.merge(
            left=pls,
            right=emb,
            how='left',
            left_index=True,
            right_index=True,
            sort=False,
        )
        logger.info('Added embedding')

        pls = pandas.merge(
            left=pls,
            right=pf,
            how='left',
            left_index=True,
            right_index=True,
            sort=False,
        )
        logger.info('Added PlasmidFinder hits')

        pls = pandas.merge(
            left=pls,
            right=pmlst,
            how='left',
            left_index=True,
            right_index=True,
            sort=False,
        )
        logger.info('Added pMLST hits')
        
        pls = pandas.merge(
            left=pls,
            right=mob,
            how='left',
            left_index=True,
            right_index=True,
            sort=False,
        )
        logger.info('Added mob hits')

        # save
        pls.to_csv(
            path_or_buf=output[0],
            sep='\t',
            header=True,
            index=False,
            index_label=False
        )




######################################################################
# GENERAL RULES
######################################################################


# BLAST DBs from FASTA
rule blastndb:
    input:
        fasta="{basename}.fna",
        bin=BIN_MAKEBLASTDB
    output:
        ["{basename}.fna.%s" % ext for ext in config['blast']['ext']]
    log:
        "{basename}.makeblastdb.log"
    params:
        date=timestamp
    message:
        "Create BLASTDB from {input.fasta}"
    shell:
        "{input.bin} -in {input.fasta} -input_type fasta -dbtype nucl -title \"plsdb_{params.date}\" -logfile {log}"

# Mash sketch
rule mash_sketch:
    input:
        fasta="{basename}.fna",
        bin=BIN_MASH
    output:
        "{basename}.msh"
    params:
        params=config['mash']['sketch_params']
    message:
        "Create Mash signatures from {input}"
    shell:
        "{input.bin} sketch {params.params} -o $(dirname {output})/$(basename -s .msh {output}) {input.fasta}"

# Mash dist (only if dist = 0)
rule mash_dist_zero:
    input:
        "{basename}.msh"
    output:
        "{basename}.dist0"
    params:
        params=config['mash']['dist0_params'],
        mash=BIN_MASH
    message:
        "Compare signatures in {input}"
    shell:
        "{params.mash} dist {params.params} {input} {input} > {output}"

# Mash dist for highly sim. seq.s (dist cutoff)
rule mash_dist_sim:
    input:
        "{basename}.msh"
    output:
        "{basename}.distS"
    params:
        params=config['mash']['distS_params'],
        mash=BIN_MASH
    message:
        "Compare signatures in {input}"
    shell:
        "{params.mash} dist {params.params} {input} {input} > {output}"

# Mash dist (all, tab format)
rule mash_dist:
    input:
        "{basename}.msh"
    output:
        "{basename}.dist"
    params:
        params=config['mash']['dist_params'],
        mash=BIN_MASH
    message:
        "Compare signatures in {input}"
    shell:
        "{params.mash} dist {params.params} {input} {input} > {output}"

# Embedding using UMAP on Mash distances
rule umap:
    input:
        "{basename}.dist"
    output:
        "{basename}.umap"
    log:
        "{basename}.umap.log"
    params:
        neighbors=config['umap']['neighbors'],
        components=config['umap']['components'],
        min_dist=config['umap']['min_dist']
    message:
        "UMAP on {input}"
    run:
        import umap, numpy

        logger = setup_logger(logging.INFO)

        # IDs
        ids = None
        with open(input[0], 'r') as ifile:
            ids = ifile.readline().rstrip('\n').split('\t')[1:]
            logger.info('Distance matrix: {num} x {num}'.format(num=len(ids)))

        # distance matrix
        logger.info('Start loading distances...')
        dist = numpy.loadtxt(
            fname=input[0],
            skiprows=1,
            usecols=range(1,len(ids) + 1) # skip 1st (contains IDs)
        )

        # embedding
        logger.info('Start embedding...')
        embedding = umap.UMAP(
            n_neighbors=params.neighbors,
            n_components=params.components,
            min_dist=params.min_dist,
            init='random',
            metric='precomputed',
            random_state=42
        ).fit_transform(dist)
        logger.info('Done.')

        # save to file
        embedding = pandas.DataFrame(
            embedding,
            columns=['D1', 'D2'],
            index=ids
        )
        embedding.to_csv(
            path_or_buf=output[0],
            sep='\t',
            header=True,
            index=True,
            index_label='ID'
        )

rule infer_host:
    input:
        #PLASMIDS_DIS_ONT
        PLASMIDS_FILT4
    output:
        PLASMIDS_HOST
    conda:
        "../envs/r_host.yml"
    log:
        os.path.join(config["logs"]["odir"]["reports"], "%s.host.log" % timestamp)
    message:
        "Inferring host for {input}"
    params:
        host_maps=config['data']['host']
    shell:
        """
        Rscript scripts/infer_host.R {input} {params.host_maps} {output}
        """

rule apply_disease_ontology:
    input:
        pls=PLASMIDS_HOST,
        dis=DISEASE_ONT,
    output:
        PLASMIDS_DIS_ONT
    conda:
        "../envs/py_env.yml"
    log:
        os.path.join(config["logs"]["odir"]["reports"], "%s.ont.log" % timestamp)
    message:
        "Applying disease ontology on {input.pls}"
    params:
        cores=96
    shell:
#        """
#        cp {input.pls} {output}
#        """
        """
        python scripts/apply_disease_ontology_concurrent.py -i {input.pls} -d {input.dis} -c {params.cores} -o {output} > {log}
        """
