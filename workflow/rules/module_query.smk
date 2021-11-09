#Module aggregating all the rules requiring an internet connection for downloading


# Plasmid query
##################################################

rule query_plasmids:
    output:
        "{dir}/{date}__{db}__nuccore.txt"
    message:
        "Query for plasmids in {wildcards.db}"
    params:
        query=lambda wildcards: config['eutils']['query']['plasmid']['cmd'].format(
            esearch_query=config['eutils']['query']['plasmid']['esearch_query'],
            efilter_query=config['eutils']['query']['plasmid']['efilter_query'].format(db=wildcards.db),
            xtract_query=config['eutils']['query']['plasmid']['xtract_query']
        ),
        header='\t'.join(config['eutils']['header']['plasmid']),
        minlength=config['filtering']['minlength']
    conda:
        "../envs/eutils.yml"
    shell:
        """
        mkdir -p $(dirname {output}) && echo -e \"{params.header}\" > {output} && {params.query} | awk -F "\\t" '{{if ($NF>{params.minlength}) print $0}}'>> {output}
        """
#Optional move to script
rule queried_plasmids:
    input:
        PLASMIDS
    output:
        PLASMIDS_ALL
    message:
        "Creating a table of all found plasmids"
    log:
        os.path.join(config["logs"]["odir"]["reports"],'%s__nuccore.txt' % timestamp)
    run:
        logger = setup_logger(logging.INFO,log[0])
        dfs = []
        for df_file in input:
            df = pandas.read_csv(df_file, sep='\t', header=0)
            df_source = df_file.split('__')[1]
            assert df_source in dbs.keys(), 'Unknown source \"{}\"'.format(df_source)
            df['Source_NUCCORE'] = df_source
            dfs.append(df)
        dfs = pandas.concat(dfs)
        dfs.to_csv(output[0], sep='\t', index=False, header=True)

# Plasmid info query
##################################################
rule query_linked_assemblies:
    input:
        PLASMIDS_ALL
    output:
        PLASMIDS_ALL_ASM
    message:
        "Query for linked Assemblies from {input}"
    params:
        header='\t'.join(config['eutils']['header']['link_asm']),
        cmd= lambda wildcards:config['eutils']['query']['linked'].replace('"','\\"')
    log:
        os.path.join(config["logs"]["odir"]["reports"],'%s__nuccore_assembly.log' % timestamp)
    shell:
         """python3 scripts/query_wrapper.py \
        --log_file {log} \
        --df_file {input} \
        --ofile {output} \
        --header \"{params.header}\" \
        --cmd \"{params.cmd}\" \
        --df_col UID_NUCCORE  \
        --split_size 200  \
        --db_source nuccore \
        --db_target assembly """

rule query_assemblies:
    input:
        PLASMIDS_ALL_ASM
    output:
        PLASMIDS_ASM
    message:
        "Query for Assemblies from {input}"
    params:
        header='\t'.join(config['eutils']['header']['assembly']),
        cmd= lambda wildcards:config['eutils']['query']['assembly'].replace('"','\\"')
    log:
        os.path.join(config["logs"]["odir"]["reports"],'%s__assembly.log' % timestamp)
    shell:
        """python3 scripts/query_wrapper.py \
        --log_file {log} \
        --df_file {input} \
        --ofile {output} \
        --header \"{params.header}\" \
        --cmd \"{params.cmd}\" \
        --df_col UID_ASSEMBLY  \
        --split_size 200  \
        --split_str ';' """

# Linked biosamples
rule query_linked_biosamples:
    input:
        PLASMIDS_ALL
    output:
        PLASMIDS_ALL_BIOS
    message:
        "Query for linked BioSamples from {input}"
    params:
        header='\t'.join(config['eutils']['header']['link_bios']),
        cmd= lambda wildcards:config['eutils']['query']['linked'].replace('"','\\"')
    log:
        os.path.join(config["logs"]["odir"]["reports"],'%s__nuccore_biosample.log' % timestamp)
    shell:
      """python3 scripts/query_wrapper.py \
      --log_file {log} \
      --df_file {input} \
      --ofile {output} \
      --header \"{params.header}\" \
      --cmd \"{params.cmd}\" \
      --df_col UID_NUCCORE  \
      --split_size 200  \
      --db_source nuccore \
      --db_target biosample"""


rule query_biosamples:
    input:
        PLASMIDS_ALL_BIOS
    output:
        PLASMIDS_BIOS
    message:
        "Query for BioSamples from {input}"
    params:
        header='\t'.join(config['eutils']['header']['biosample']),
        cmd= lambda wildcards:config['eutils']['query']['biosample'].replace('"','\\"')
    log:
        os.path.join(config["logs"]["odir"]["reports"],'%s__biosample.log' % timestamp)
    shell:
        """python3 scripts/query_wrapper.py \
        --log_file {log} \
        --df_file {input} \
        --ofile {output} \
        --header \"{params.header}\" \
        --cmd \"{params.cmd}\" \
        --df_col UID_BIOSAMPLE  \
        --split_size 300  \
        --split_str ';' && \
        sed -i  's/\\t[nN]ot [cC]ollected/\\tNA/g' {output} &&\
        sed -i  's/\\t[nN]ot [aA]pplicable/\\tNA/g' {output} &&\
        sed -i  's/\\t[nN]ot [dD]etermined/\\tNA/g' {output} &&\
        sed -i  's/\\t[nN]ot [aA]vailable/\\tNA/g' {output} &&\
        sed -i  's/\\t[mM]issing/\\tNA/g' {output} &&\
        sed -i  's/\\t[uU]nknown/\\tNA/g' {output}
        """
#A better solution than this sed chain would be desirable

#TODO move to script
# Taxon query
##################################################
rule query_taxon:
    input:
        PLASMIDS_ALL
    output:
        PLASMIDS_TAX
    message:
        "Query for linked Taxa from {input}"
    run:
        logger = setup_logger(logging.INFO)

        ids = pandas.read_csv(input[0], sep='\t', header=0, dtype=str)['TaxonID_NUCCORE'].fillna(value='')
        ids = set(ids.values)
        assert '' not in ids, 'Have an empty ID in {}'.format(input[0])
        logger.info('There are {} unique taxonomy IDs'.format(len(ids)))

        cmd = "{query} > {ofile}".format(
            query=config['eutils']['query']['taxon'].format(
                ids=','.join(ids)
            ),
            ofile=output[0]
        )
        cmd, cmd_s, cmd_o = run_cmd(cmd)
        assert cmd_s == 0, 'CMD: {}: {}\n{}'.format(cmd, cmd_s, cmd_o)

#TODO move to script
# Location query
##################################################
rule parse_locations:
    input:
        PLASMIDS_FILT1
    output:
        PLASMIDS_FILT1_LOC
    message:
        "Parsing locations from {input}"
    params:
        api_keys=config['data']['api_keys'],
        locs=config['data']['locs']
    run:
        from utils import preproc_loc_str, preproc_loc_coords, parse_location
        from utils import load_locs, save_locs, update_locs
        logger = setup_logger(logging.INFO)

        # api key
        with open(params.api_keys, 'r') as ifile:
            api_key = [k.rstrip('\n').strip() for k in ifile.readlines()]
            api_key = api_key.pop()
        logger.info('Using API key: {}'.format(api_key))

        # load known locations
        locs = load_locs(params.locs)

        # sample table with locations
        pls = pandas.read_csv(input[0], header=0, sep='\t', dtype=str)
        pls.set_index('UID_NUCCORE', drop=False, inplace=True, verify_integrity=True)

        # parse
        coords = []
        for i in tqdm(pls.index):
            # location name and coordinates
            l_n = pls.loc[i,'Location_BIOSAMPLE']
            l_c = pls.loc[i,'Coordinates_BIOSAMPLE']

            # no location data
            if l_n is None and l_c is None:
                continue

            # pre-processing
            l_n = preproc_loc_str(l_n)
            l_c = preproc_loc_coords(preproc_loc_str(l_c))
            if l_n is None and l_c is None:
                continue

            # at least name or coordinates given
            if l_c is not None: # location coordinates
                l_c_str = '{};{}'.format(l_c[0], l_c[1])
                if locs is None or l_c_str not in locs['location']:
                    logger.info('Retrieving coordinates for location: \"{}\"'.format(l_c))
                    try:
                        parsed = parse_location(loc_str=l_c, api_key=api_key, is_name=False)
                        locs = update_locs(locs, {'location': l_c_str, 'type': 'coordinates', 'lat': parsed['lat'], 'lng': parsed['lng']})
                        save_locs(locs, params.locs)
                    except Exception as e:
                        logger.error('Error while retrieving coordinates for location \"{}\"'.format(l_c))
                        raise(e)
                coords.append({'ID': i, 'loc_lat': locs.loc[l_c_str,'lat'], 'loc_lng': locs.loc[l_c_str,'lng']})
            elif l_n is not None: # location name
                if locs is None or l_n not in locs['location']:
                    logger.info('Retrieving coordinates for location: \"{}\"'.format(l_n))
                    try:
                        parsed = parse_location(loc_str=l_n, api_key=api_key, is_name=True)
                        locs = update_locs(locs, {'location': l_n, 'type': 'name', 'lat': parsed['lat'], 'lng': parsed['lng']})
                        save_locs(locs, params.locs)
                    except Exception as e:
                        logger.error('Error while retrieving coordinates for location \"{}\"'.format(l_n))
                        raise(e)
                coords.append({'ID': i, 'loc_lat': locs.loc[l_n,'lat'], 'loc_lng': locs.loc[l_n,'lng'], 'loc_parsed': l_n})
        coords = pandas.DataFrame(coords)
        coords.set_index('ID', drop=True, verify_integrity=True, inplace=True)
        pls = pandas.merge(
            left=pls,
            right=coords,
            how='left',
            left_index=True,
            right_index=True,
            sort=False,
        )

        # save
        pls.to_csv(output[0], sep='\t', header=True, index=False, index_label=False)


#TODO move to script
rule get_pmlst_data:
    input:
        PMLST_CLEAN
    output:
        PMLST_ALLELES
    message:
        "Download pMLST allele FASTA files"
    params:
        db=config['pmlst']['db'],
        url=config['pmlst']['url'],
        url_schemes=config['pmlst']['url_schemes'],
    run:
        # reference: https://github.com/kjolley/BIGSdb.git: scripts/rest_examples/python/download_alleles.py
        import requests
        from glob import glob
        from scripts.utils import proc_mlst_scheme_name, download_pmlst_scheme_alleles, download_pmlst_scheme_profiles
        
        logger = setup_logger(logging.INFO)

        # get schemes
        schemes = requests.get(params.url_schemes)
        assert schemes.status_code != 404, 'Invalid URL {}'.format(params.url_schemes)
        schemes = schemes.json()['schemes']
        logger.info('There are {} pMLST schemes'.format(len(schemes)))

        # get loci
        for scheme in schemes:
            # scheme name and URL
            scheme_name = scheme['description']
            scheme_url  = scheme['scheme']
            scheme_name2 = proc_mlst_scheme_name(scheme_name)
            logger.info('Scheme {} ({}): {}'.format(scheme_name, scheme_name2, scheme_url))

            # where to save files
            scheme_dir = os.path.join(os.path.dirname(output[0]), scheme_name2)
            mkdir(scheme_dir, True)

            # where to save profiles
            scheme_profiles = os.path.join(scheme_dir, scheme_name2 + '.txt')

            # alleles
            download_pmlst_scheme_alleles(scheme_name, scheme_url, scheme_dir)

            # profiles
            profiles_url = params.url +  '/db/' + params.db + '/schemes/' + scheme_url.split('/')[-1]
            download_pmlst_scheme_profiles(scheme_name, profiles_url, scheme_dir, scheme_profiles)

        with open(output[0], 'w') as ofile:
            ofile.write('done')

ruleorder: xtract_fasta2 > get_fasta
# FASTA for nuccore UIDs in table
rule get_fasta:
    input:
        "{basename}.txt"
    output:
        "{basename}.fna"
    message:
        "Downloading sequences for IDs in {input}"
    shell:
        """
        mkdir -p $(dirname {output}) &&
        python scripts/download_fastas.py -t {input} -i \"UID_NUCCORE\" -o {output} -c 1 -s 500 -e $(dirname $(which efetch)) 2>&1 | tee {output}.log &&
        cat {output}.tmp.* > {output} && rm {output}.tmp.*
        """

rule xtract_fasta2:
    input:
        pls=PLASMIDS_FILT2,
        fna=PLASMIDS_FILT1_FASTA
    output:
        PLASMIDS_FILT2_FASTA 
    message:
        "Extracting sequences from {input.fna}"
    log:
        os.path.join(config['logs']['odir']['reports'], '%s__filtered2.fna.log' % timestamp)
    shell:
        """
        python scripts/subsample_fasta.py -i {input.pls} -f {input.fna} -o {output} -l {log}
        """

rule xtract_fasta3:
    input:
        pls=PLASMIDS_FILT3,
        fna=PLASMIDS_FILT2_FASTA
    output:
        PLASMIDS_FILT3_FASTA 
    message:
        "Extracting sequences from {input.fna}"
    log:
        os.path.join(config['logs']['odir']['reports'], '%s__filtered4.fna.log' % timestamp)
    shell:
        """
        python scripts/subsample_fasta.py -i {input.pls} -f {input.fna} -o {output} -l {log}
        """

rule download_dis_ont:
    output:
        DISEASE_ONT
    message:
        "Downloading disease ontology"
    shell:
        """
        wget http://purl.obolibrary.org/obo/doid.obo  -O {output}
        """
