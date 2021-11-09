
# Plasmid filtering (1)
##################################################
rule filter1:
    input:
        PLASMIDS_FULL1
    output:
        PLASMIDS_FILT1
    message:
        "Filter plasmids from {input}"
    params:
        dfilter=config['filtering']['dfilter']
    log:
        os.path.join(config['logs']['odir']['reports'], '%s__filtered1.log' % timestamp)
    run:
        logger = setup_logger(logging.INFO,log[0])
        pls = pandas.read_csv(input[0], sep='\t', header=0, dtype=str)
        pls.set_index('UID_NUCCORE', drop=False, inplace=True, verify_integrity=True)
        logger.info('Read in {} plasmid records'.format(pls.shape[0]))

        # filter by description
        keep_rows = pandas.Series([re.search(params.dfilter, d, flags=re.IGNORECASE) is None for d in pls['Description_NUCCORE']], index=pls.index)
        pls = pls.loc[keep_rows,:]
        logger.info('Filtered: description "{}": kept {} records'.format(params.dfilter, pls.shape[0]))

        # filter by (assembly) completeness
        def pls_filter(uid):

            assembly_tag = pls.loc[uid, 'Status_ASSEMBLY']

            has_assembly_tag = pandas.notnull(pls.loc[uid, 'UID_ASSEMBLY']) and pandas.notnull(assembly_tag) and assembly_tag != ""

            complete_tag  = pls.loc[uid, 'Completeness_NUCCORE']

            has_complete_tag = pandas.notnull(complete_tag) and complete_tag != ""

            if has_assembly_tag and has_complete_tag:
                return assembly_tag == 'Complete Genome' and complete_tag == 'complete'
            elif not has_assembly_tag:
                return complete_tag == 'complete'
            elif not has_complete_tag:
                return assembly_tag == 'Complete Genome'
            return False

        logger.info('Completeness nuccore tag: {}'.format(set(pls['Completeness_NUCCORE'])))
        logger.info('Completeness assembly tag: {}'.format(set(pls['Status_ASSEMBLY'])))
        pls = pls.loc[pls.index.map(pls_filter),:]
        logger.info('Filtered: (assembly) completeness: kept {} records'.format(pls.shape[0]))

        # filter by taxonomy
        pls = pls.loc[pls['taxon_superkingdom_id'] == '2',:]
        logger.info('Filtered: superkingdom ID is 2 (i.e. Bacteria): kept {} records'.format(pls.shape[0]))

        # save
        pls.to_csv(output[0], sep='\t', header=True, index=False)

# Plasmid filtering (2)
##################################################

rule filter2:
    input:
        pls=PLASMIDS_FILT1_LOC,
        fna=PLASMIDS_FILT1_FASTA,
        dist=PLASMIDS_FILT1_DIST0
    output:
        PLASMIDS_FILT2
    message:
        "Filter out identical plasmids {input.pls}"
    log:
        os.path.join(config['logs']['odir']['reports'], '%s__filtered2.log' % timestamp)
    run:
        import numpy
        from Bio import SeqIO
        from Bio.SeqUtils import GC
        from scripts.utils import str2timestamp

        logger = setup_logger(logging.INFO,log[0])

        pls = pandas.read_csv(input.pls, sep='\t', header=0, dtype=str)
        pls.set_index('ACC_NUCCORE', drop=False, inplace=True, verify_integrity=True)
        logger.info('Read in {} plasmid records'.format(pls.shape[0]))

        # collect seq. stats.s
        logger.info('Collecting sequences and their data ...')
        seq = []
        seqs = {}
        with open(input.fna, 'r') as ifile:
            for record in tqdm(SeqIO.parse(ifile, 'fasta')):
                assert re.fullmatch(r'\w+\.\d+', record.id), 'Unexpected ID format in FASTA: {}'.format(record.id)
                assert record.id in list(pls.index), 'Unknown FASTA ID {}'.format(record.id)
                seq.append({
                    'ACC_NUCCORE': record.id,
                    'GC_NUCCORE': GC(record.seq),
                    'Length': len(record.seq)
                })
                seqs[record.id] = record.seq
        seq = pandas.DataFrame(seq)
        seq.set_index('ACC_NUCCORE', drop=True, inplace=True, verify_integrity=True)
        # add to table
        pls = pandas.merge(
            left=pls,
            right=seq,
            how='left',
            left_index=True,
            right_index=True,
            sort=False,
        )
        # check seq. length
        assert pls['Length_NUCCORE'].astype(int).equals(pls['Length'].astype(int)), 'Sequence length is not always identical: nuccore vs. FASTA'
        # drop 2nd seq. len. col.
        pls.drop(columns=['Length'], inplace=True)

        # find and assign identity groups
        logger.info('Grouping identical sequences ...')
        pls['UID_IDENTGROUP'] = None
        groupID = 0
        with open(input.dist, 'r') as ifile:
            for line in tqdm(ifile):
                sID, qID, dist, pv, sh = line.rstrip('\n').split('\t')
                # same ID -> skip
                if sID == qID:
                    continue
                # group ID already set -> skip
                if pandas.notnull(pls.loc[sID, 'UID_IDENTGROUP']) and pandas.notnull(pls.loc[qID, 'UID_IDENTGROUP']):
                    continue
                # check if really equal
                if len(seqs[sID]) != len(seqs[qID]):
                    continue
                elif seqs[sID] != seqs[qID]: # length is equal -> check sequences
                    continue
                # sequences are identical
                # one of the seq.s has already an ID -> set to same value
                if pandas.notnull(pls.loc[sID, 'UID_IDENTGROUP']):
                    pls.loc[qID, 'UID_IDENTGROUP'] = pls.loc[sID, 'UID_IDENTGROUP']
                    continue
                if pandas.notnull(pls.loc[qID, 'UID_IDENTGROUP']):
                    pls.loc[sID, 'UID_IDENTGROUP'] = pls.loc[qID, 'UID_IDENTGROUP']
                    continue
                # new group -> set both to new value
                groupID += 1
                pls.loc[sID, 'UID_IDENTGROUP'] = groupID
                pls.loc[qID, 'UID_IDENTGROUP'] = groupID
        # save
        pls.to_csv(output[0] + '.backup', sep='\t', header=True, index=False, index_label=False)

        # filter
        def compare_records(rec1, rec2):
            """
            Compare two plasmid records:
            0) Prefer the one with higher version number (if same accession)
            1) Prefer the one from RefSeq
            2) Prefer the one with location information
            3) Prefer the one with an assembly
            4) Prefer the older one
            5) Compare by accession (string), prefer the "smaller" one (just have something to break the ties)
            """
            rec1_acc, rec1_ver = rec1['ACC_NUCCORE'].split('.')
            rec1_ver = int(rec1_ver)
            rec2_acc, rec2_ver = rec2['ACC_NUCCORE'].split('.')
            rec2_ver = int(rec2_ver)

            if rec1_acc == rec2_acc:
                assert rec1_ver != rec2_ver, 'Both records have the same accession and version: {} vs {}'.format(rec1['UID_NUCCORE'], rec1['UID_NUCCORE'])
                if rec1_ver > rec2_ver:
                    return rec1
                else:
                    return rec2
            elif rec1['Source_NUCCORE'] == "refseq" and rec2['Source_NUCCORE'] != "refseq":
                return rec1
            elif rec1['Source_NUCCORE'] != "refseq" and rec2['Source_NUCCORE'] == "refseq":
                return rec2
            elif pandas.notnull(rec1['loc_lat']) and pandas.isnull(rec2['loc_lat']):
                return rec1
            elif pandas.isnull(rec1['loc_lat'])  and pandas.notnull(rec2['loc_lat']):
                return rec2
            elif pandas.notnull(rec1['UID_ASSEMBLY']) and pandas.isnull(rec2['UID_ASSEMBLY']):
                return rec1
            elif pandas.isnull(rec1['UID_ASSEMBLY'])  and pandas.notnull(rec2['UID_ASSEMBLY']):
                return rec2
            elif str2timestamp(rec1['CreateDate_NUCCORE'], '%Y/%m/%d') < str2timestamp(rec2['CreateDate_NUCCORE'], '%Y/%m/%d'):
                return rec1
            elif str2timestamp(rec1['CreateDate_NUCCORE'], '%Y/%m/%d') > str2timestamp(rec2['CreateDate_NUCCORE'], '%Y/%m/%d'):
                return rec2
            else:
                if rec1['ACC_NUCCORE'] < rec2['ACC_NUCCORE']:
                    return rec1
                else:
                    return rec2

        logger.info('Filter by identity ...')
        pls['Identical'] = None
        keep = pandas.Series(True, index=pls['ACC_NUCCORE'])
        for gr, gr_df in tqdm(pls.loc[pandas.notnull(pls['UID_IDENTGROUP']),:].groupby(by=['UID_IDENTGROUP'])):
            assert pandas.notnull(gr)
            best = None
            for acc in gr_df.index:
                if best is None:
                    best = acc
                else:
                    best = compare_records(gr_df.loc[best,:], gr_df.loc[acc,:])['ACC_NUCCORE']
            keep.loc[gr_df.loc[gr_df['ACC_NUCCORE'] != best,'ACC_NUCCORE']] = False
            pls.loc[best,'Identical'] = ';'.join(sorted([acc for acc in gr_df['ACC_NUCCORE'] if acc != best]))
        pls = pls.loc[keep.index[keep]]
        pls.drop(columns=['UID_IDENTGROUP'], inplace=True)
        logger.info('Filtered: indentical sequences: kept {} records'.format(pls.shape[0]))

        # by accession and version
        logger.info('Filter by accession ...')
        pls['ACC_WO_VERSION'] = pls['ACC_NUCCORE'].apply(lambda x: x.split('.')[0])
        pls['ACC_WO_VERSION'] = pls['ACC_WO_VERSION'].apply(lambda x: x.split('_')[1] if '_' in x else x)
        pls['OldVersion'] = None
        keep = pandas.Series(True, index=pls['ACC_NUCCORE'])
        for gr, gr_df in tqdm(pls.groupby(by=['ACC_WO_VERSION'])):
            assert pandas.notnull(gr)
            best = None
            for acc in gr_df.index:
                if best is None:
                    best = acc
                else:
                    best = compare_records(gr_df.loc[best,:], gr_df.loc[acc,:])['ACC_NUCCORE']
            keep.loc[gr_df.loc[gr_df['ACC_NUCCORE'] != best,'ACC_NUCCORE']] = False
            pls.loc[best,'OldVersion'] = ';'.join(sorted([acc for acc in gr_df['ACC_NUCCORE'] if acc != best]))
        pls = pls.loc[keep.index[keep]]
        pls.drop(columns=['ACC_WO_VERSION'], inplace=True)
        logger.info('Filtered: sequences w/ same accession: kept {} records'.format(pls.shape[0]))

        # save
        pls.to_csv(output[0], sep='\t', header=True, index=False, index_label=False)

# Plasmid filtering (3)
##################################################
rule filter3:
    input:
        pls=PLASMIDS_FILT2,
        fna=PLASMIDS_FILT2_FASTA,
        rmlst=PLASMIDS_FILT2_RMLST
    output:
        PLASMIDS_FILT3
    message:
        "Filter based on rMLST results"
    params:
        cores=30,
        cutoff=config['rmlst']['rmlst_max_loci'],
        blastn_bin=BIN_BLASTN,
        blastn_header=config['rmlst']['blastn_header'],
        blastn_pident=config['rmlst']['blastn_pident'],
        blastn_qcovs=config['rmlst']['blastn_qcovs']
    log:
        os.path.join(config['logs']['odir']['reports'], '%s__filtered3.log' % timestamp)
    run:
        from multiprocessing import Pool

        logger = setup_logger(logging.INFO,log[0])

        pls = pandas.read_csv(input.pls, sep='\t', header=0, dtype=str)
        pls.set_index('ACC_NUCCORE', drop=False, inplace=True, verify_integrity=True)
        logger.info('Read in {} plasmid records\n{}'.format(pls.shape[0], pls.head()))

        rmlst = pandas.read_csv(input.rmlst, sep='\t', header=0, dtype=str)
        logger.info('Read in {} rMLST hits\n{}'.format(rmlst.shape[0], rmlst.head()))

        # count unique loci per query
        pls['hits_rMLST'] = ''
        pls['hitscount_rMLST'] = 0
        for qseqid, df in rmlst[['qseqid', 'slocus']].groupby(['qseqid']):
            hits = sorted(list(set(df['slocus'])))
            pls.loc[qseqid, 'hits_rMLST'] = ';'.join(hits)
            pls.loc[qseqid, 'hitscount_rMLST'] = len(hits)
        pls.to_csv('%s.backup' % output[0], sep='\t', header=True, index=False, index_label=False)

        # check those above given cutoff
        check_ids = pls.loc[pls['hitscount_rMLST'] > params.cutoff,'ACC_NUCCORE']
        # run BLASTn
        from scripts.utils import run_blastn_check
        pool = Pool(params.cores)
        tmp = pool.starmap(
            run_blastn_check,
            [(acc, output[0], input.fna, config['rmlst']['blastn'], params.blastn_bin, params.blastn_header, params.blastn_pident) for acc in check_ids]
        )
        pool.close()
        pool.join()
        
        # check results
        for acc in check_ids:
            # clean up: rm FASTA
            acc_fasta = '%s.%s.fna' % (output[0], acc)
            if os.path.exists(acc_fasta):
                os.remove(acc_fasta)
            # load hits
            acc_ofile = '%s.%s.tsv' % (output[0], acc)
            try:
                acc_df = pandas.read_csv(acc_ofile, sep='\t', header=None, names=params.blastn_header)
            except pandas.errors.EmptyDataError:
                logger.info('KEEP {}: no blastn hits'.format(acc))
                continue
            # filter hits
            acc_df = acc_df.loc[(acc_df['pident'] >= params.blastn_pident) & (acc_df['qcovs'] >= params.blastn_qcovs),:]
            # no hits left
            if acc_df.shape[0] == 0:
                logger.info('KEEP {}: no blastn hits left after filtering (pident {}, qcovs {})'.format(acc, params.blastn_pident, params.blastn_qcovs))
                continue
            else:
                top1 = acc_df.iloc[0,:]
                logger.info('RM {}: 1st hit: {} [{}] (evalue {}, pident {}, qcovs {})'.format(acc, top1['sseqid'], top1['stitle'], top1['evalue'], top1['pident'], top1['qcovs']))
                pls.drop(index=[acc], inplace=True)
        logger.info('Filtered: by number of rMLST and blastn hits: kept {} records'.format(pls.shape[0]))

        # save
        pls.to_csv(output[0], sep='\t', header=True, index=False, index_label=False)


rule filter4:
    input:
        PLASMIDS_FILT3, 
        PLASMIDS_FILT3_FASTA
    output:
        PLASMIDS_FILT4
    message:
        "Remove artifacts from plasmid collection"
    log:
        os.path.join(config['logs']['odir']['reports'], '%s__filtered5.log' % timestamp)
    shell:
        """
        python scripts/filt4_remove_artifacts.py -i {input[0]} -f {input[1]} -l {log[0]} -o {output[0]}
        """

