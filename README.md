
Pipeline for data collection
==============================

## News
Our manuscript discussing the new features of PLSDB was accepted to the annual 2022 Nucleic Acid Research database Issue! The manuscript can be found [here](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab1111/6439675).

## Summary
![pipeline graph](workflow/dag.svg)

- `retrieve_plasmid`: Retrieve plasmid data from NCBI database
    - - `esearch` query from **Orlek et al.** 
    - Used databases: nucleotide, biosample, assembly, and taxonomy database. NCBI database sources:
        - INSDC (DDBJ, EMBL/ENA, GenBank): [International Nucleotide Sequence Database Collaboration](https://www.ncbi.nlm.nih.gov/genbank/collab/)
        - RefSeq
- `filter_metadata`: To remove incomplete or nonbacterial records these are filtered by metadata attributes:
    - Duplicated entry (NUCCORE_DuplicatedEntry)
    - Record description (regular expression from **Orlek et al.**)
    - Assembly Lastest (True)
    - Assembly Completeness (Complete)
        - If no assembly: Completeness status of the nuccore record has to be `complete`
        - Has assembly: assembly status of the latest version has to be `Complete genome`
            - [NCBI assembly help page](https://www.ncbi.nlm.nih.gov/assembly/help/)
    - By taxonomy: superkingdom taxon should be Bacteria
- `filter_sequences`: To remove identical records
    - Group plasmids with identical sequences
    - Among these groups select one record
        - Prefer the one from RefSeq
        - Prefer the one with location information
        - Prefer the one with an assembly
        - Prefer the most recent assembly release date
        - Prefer the one with a biosample
        - Prefer the newewst nuccore creation date
        - Prefer the one with the highest coverage
        - If all equals, choose the first one
- `filter_rmlst`: To remove putative chromosomal sequences
    - Obtain rMLST database
    - Create a local NCBI chromosomal sequences database
    - The plasmid sequences are aligned against the rMLST allele sequences and local NCBI-db
    - Records having more than 5 unique rMLST loci are searched in NCBI chromosomal sequences using BLASTn (remote access)
    - Records with hits are removed
- `filter_artifacts`: Remove possible artifacts sequences
- `process_abricate`: Annotate antimicrobial resistance or virulence genes.
    - BLASTn search in DBs provided by ABRicate
        - Blaster from [CGE core module](https://bitbucket.org/genomicepidemiology/cge_core_module) is used for search and pre-processing
        - Filtering:
            - Identity and coverage cutoffs
            - Overlapping matches are removed
        - All hits are ollected into one file
- `process_pmslt`: Annotate using pMLST
    - For each found replicon use the associated pMLST scheme (if available)
    - Use `mlst` to perform the pMLST analysis
    - Process the results
        - Set IncF ST according to the FAB formula (**Villa et al.**)
    - Create BLAST database file from plasmid FASTA
    - Create sketches from plasmid FASTA using Mash
- `dstream_sim_records`: List of similar plasmids
    - Use Mash do compute pairwise distances (use a distance cutoff)
    - Create a list of unique pairs
- `dstream_umap`: Embedding
    - Compute pairwise distances between plasmids using Mash
    - Compute embedding using UMAP
- `process_infotable`: Create info table
    - Record information
    - Embedding coordinates
    - PlasmidFinder hits
    - pMLST hits
- `dstream_compare`: Compare created table to an older version
    - Which plasmid records were removed
    - Which plasmid records were added
    - Which plasmid records changed

## Preparations

### PubMLST data
This data processing pipeline makes use of the [PubMLST website](https://pubmlst.org/) developed by Keith Jolley ([Jolley & Maiden 2010, BMC Bioinformatics, 11:595](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-595)) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust.
#### rMLST data
*Note: requires a PubMLST account*
*Note: requires graphical interface*
*Note: PubMLST account needs to request access to Ribosomal MLST locus/sequence definitions database from rMLST admin. Access is normally granted within a day. *
*Error: `Message: 'chromedriver' executable needs to be in PATH` : Make sure that chromedriver is installed. 
As `chromium` already comes with a `chromedriver` installation, you can try: `sudo pacman -Syyu` followed by `sudo pacman -S chromium`*

To remove putative chromosomal sequences rMLST analysis is performed which requires rMLST sequences from [PubMLST](https://pubmlst.org/rmlst/). There is an API for the PubMLST services, however using it seems to require much more effort than downloading the data through a web browser. Thus, there is a rule ([retrieve_rmlst_data](workflow/rules/module_retrival.smk)) that downloads the sequences automatically (given the login data). **This rule needs a graphical interface**, please, run this rule locally in your computer.
  
Here, a login and password are required. **Please, create and account and specify your credentials in `config.yml`**.

*Note: Cookie agreement might cause problems. Requires minor changes if "Got it!" is changed to different link text.*

#### pMLST
There is a mapping from PlasmidFinder IDs to pMLST profile names in `pipeline.json` (`pmlst/map`).
It may require an update depending on which pMLST schemes are available from PubMLST and which IDs are currently in the PlasmidFinder database.
- *Note: Information on pMLST schemes is shown in: https://pubmlst.org/plasmid/*

- Path to installed pMLST schemes: `~/miniconda3/envs/plsdb/db/pubmlst/`
    - Each scheme is one directory (also listed in the log file when created)
- Path to installed PlasmidFinder DB: `~/miniconda3/envs/plsdb/db/plasmidfinder/sequences`

### ABRicate
Please, if the most recent version of ABRicate contains the most recent database links [abricate-get_db](https://github.com/tseemann/abricate/blob/master/bin/abricate-get_db).
- [VFDB](http://www.mgc.ac.cn/VFs/download.htm)
- [CARD](https://card.mcmaster.ca/download)
- [ARG-ANNOT](https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/)
- [Plasmidfinder](https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src/master/)
- [Resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/)

**IMPORTANT**: Currently, ABRicate (version `1.0.1`) does not update some databases correctly:
- ARG-ANNOT: The URL changed
The file `patch_abricate-get_db` (`abricate: getdb_bin` param in [config.yml](workflow/config.yml)) should resolve these issues. It is a copy of the `abricate-get_db`, but changing the URL of ARG-ANNOT.  If you also find some deprecated links, please substitute them and set the param `abricate: replace_getdb: True` in [config.yml](workflow/config.yml).

### API keys
#### NCBI data
To retrieve data from NCBI, please obtain an [API](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us) and specify it in the [`config.yml`](workflow/config.yml).
#### Location queries
To map location names to coordinates the Nominatim API and Google API are used. Google API is only used for comparative purposes, as their policy doesn't allow the 
storage of google's content (more [here](https://developers.google.com/maps/documentation/geocoding/policies)).Google requires API key, please which requires you register ([more info](https://support.google.com/googleapi/answer/6158862?hl=en)).

### BIOSAMPLE_Host
Already known hosts are in [hosts_version.csv](src/hosts_20231101.csv).
Run the rule [`process_create_host_mapping`](workflow/rules/module_process.smk) and manually check the new versions of `host` mapping. Find more details at the end of the log file of the rule (`logs/process_create_host_mapping.log`) or in the rule [`process_manually_inspect_hosts`](workflow/rules/module_process.smk).

### BIOSAMPLE_Location
Already known locations are in [locations_version.csv](src/locations_20231101_v2.csv) and corrections to find some specific locations are display in [location_correction_version.csv](src/location_corrections_20231101_v2.csv)

Run the rule [`process_parse_locations`](workflow/rules/module_process.smk) and manually check the new versions of `location` and `location_corrections`. Find more details in the rule [`process_manually_inspect_locations`](workflow/rules/module_process.smk).


### Conda & Snakemake
If needed, install (mini-)conda
```bash
cd ~
# get miniconda (for linux)
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# install
bash Miniconda3-latest-Linux-x86_64.sh
# set path to binaries in your ~/.bashrc
export PATH=$HOME/miniconda3/bin:$PATH
```

If needed, install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
```bash
# If required, install mamba package manager in you base env
conda install  -c conda-forge mamba
# Install snakemake
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Current versions:  
    - conda==23.7.3
    - snakemake==7.32.4


### Comparing new and old versions
The last rule in the pipeline requires a "master" table from an older version.
The path has to be set in `config.yml` (attribute `previous_table`) and the file must exist.

## Running the pipeline

- Do **not ignore** log files. They contain timestamps and useful information to detect possible errors. Add more information if necessary.
- Take notes: Write down what was changed: version updates, big fixes etc.
- Do **not** run the pipeline on the first day of a month (usually many requests fail or return nothing)
- Do **not** execute the complete pipeline: use `snakemake --use-conda -c {CORES} --force target_rule` to run groups of rules
- Before running something, **list the commands** to be executed: `snakemake -np`. See [Section](#groups-of-execution-sequentially)
- If a step where some data is fetched from NCBI **fails** try to **re-run** the step
    - Sometimes NCBI return an empty result or the request fails
    - Re-running the command usually solves the problem
- The **longest steps** are:
    - Automatic rules:
        - retrieve_nuccoredb_seqs: Download of putative chromosomes from NCBI server to create local db (2023-10-06: ~7h)
        - process_rmlst_blastn: BLASTn search against rmlst database (2023-10-05: ~16h)
        - filtering_rmlst: filter plasmid data using information form nuccoredb and rmlst_blastn (2023-10-06: ~2h)
    - Manual curation rules:
        - process_manually_inspect_hosts: Depends on the number of unknown hosts, but save at least 2 days (16h)
        - process_manually_inspect_locations: Depends on the number of unknown locations, but save at least 1 day (8h)
- Other steps require usually only a few minutes and should run under one hour
- Try to run all steps requiring updated data on the **same day**
    - I.e. getting new data for rMLST, `abricate`, pMLST, NCBI data (retrieve rules)

### Groups of execution (sequentially)
- RETRIVAL
    - Plasmid NCBI retrival: `retrieve_plasmid_metadata retrieve_plasmid_taxid filter_metadata retrieve_fasta`
    - NCBI chromosomal sequences retrival: `retrieve_nuccoredb_ids retrieve_nuccoredb_seqs process_make_nuccoredb_blastdb`
    - ABRicate retrival: `retrieve_abricate_getdb retrieve_abricatedb`
    - pmlst retrival: `retrieve_pmlst_data process_make_pmlst_blastdb`
    - rmlst retrival (graphical interface): `retrieve_rmlst_data process_make_rmlst_blastdb`
    - human disease ontology: `retrieve_disease_ont`
- FILTERING:
    - `filter_sequences process_rmlst_blastn filter_rmlst filter_artifacts`
- PROCESSING:
    - `process_calculate_GC process_abricate process_join_abricate process_pmlst`
    - `process_mash_sketch process_mash_dist process_umap process_mash_dist_sim process_dstream_sim_records`
- MANUAL_CURATION:
    - `process_create_host_mapping process_manually_inspect_hosts process_infer_host`
    - `process_disease_ont` (intermediate step, not manual curation)
    - `process_parse_locations process_manually_inspect_locations`
- DSTREAM:
    - `process_infotable dstream_krona_xml dstream_krona_html dstream_summary dstream_compare`

# References

- **Mash**: "Mash: fast genome and metagenome distance estimation using MinHash", B. D. Ondov, T. J. Treangen, P. Melste
d, A. B. Mallonee, N. H. Bergman, S. Koren and A. M. Phillippy, Genome Biology, 2016, [paper link](https://genomebiology
.biomedcentral.com/articles/10.1186/s13059-016-0997-x), [repository link](https://github.com/marbl/Mash)
- **UMAP**: "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction", L. McInnes and J. Healy, N. Saul and L. Großberger, Journal of Open Source Software
v, 2018,
[paper link](http://joss.theoj.org/papers/10.21105/joss.00861), [repository link](https://github.com/lmcinnes/umap)
- **BLAST**: "Basic local alignment search tool." , S.F. Altschul, W. Gish, W. Miller, E. W. Myers and D. J.  Lipman, J. Mol. Biol. 215:403-410, [BLAST paper link](https://www.ncbi.nlm.nih.gov/pubmed/2231712), [BLAST+ paper link](https://www.ncbi.nlm.nih.gov/pubmed/20003500), [tool link](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+)
- **ABRicate**: Tool implemented by Thorsten Seemann [repository link](https://github.com/tseemann/abricate)
- **ARG-ANNOT**: "ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes", S. K. Gupta, B. R. Padmanabhan, S. M. Diene, R. Lopez-Rojas,
M. Kempf, L. Landraud, and J. M. Rolain, Antimicrob. Agents Chemother., 2014, [paper link](https://www.ncbi.nlm.nih.gov/pubmed/24145532)
- **CARD**: "CARD 2017: expansion and model-centric curation
of the comprehensive antibiotic resistance database.", B. Jia, A. R. Raphenya, B. Alcock, N. Waglechner, P. Guo, K. K. Tsang, B. A. Lago, B. M. Dave, S. Pereira, A. N. Sharma, S. Doshi, M. Courtot, R. Lo, L. E. Williams, J. G. Frye, T. Elsayegh, D. Sardar, E. L. Westman, A. C. Pawlowski, T. A. Johnson, F. S. Brinkman, G. D. Wright, and A. G. McArthur, Nucleic Acids Res., 2017, [paper link]( https://www.ncbi.nlm.nih.gov/pubmed/27789705)
- **ResFinder**: "Identification of acquired antimicrobial resistance genes", E. Zankari, H. Hasman, S. Cosentino, M. Vestergaard, S. Rasmussen, O. Lund, F. M. Aarestrup, and M. V. Larsen, J. Antimicrob. Chemother., 2012, [paper link](https://www.ncbi.nlm.nih.gov/pubmed/22782487)
- **VFDB**: "VFDB: a reference database for bacterial virulence factors", L. Chen, J. Yang, J. Yu, Z. Yao, L. Sun, Y. Shen, and Q. Jin, Nucleic Acids Res., 2005, [paper link](https://www.ncbi.nlm.nih.gov/pubmed/15608208)
- **PlasmidFinder**: "In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing.", A. Carattoli, E. Zankari, A. Garcia-Fernandez, M. Voldby Larsen, O. Lund, L. Villa, F. Møller Aarestrup, and H. Hasman, Antimicrob. Agents Chemother., 2014, [paper link](http://aac.asm.org/content/58/7/3895.long), [repository link](https://bitbucket.org/genomicepidemiology/plasmidfinder)
- **pMLST in PubMLST**: [web-site](https://pubmlst.org/plasmid/)
- **mlst**: Tool implemented by Thorsten Seemann, [repository link](https://github.com/tseemann/mlst)
- **OpenCageData**: An API to convert coordinates to and from places, [web-site](https://opencagedata.com/)
- **rMLST**: [rMLST at PubMLST](https://pubmlst.org/rmlst/)
- **Jolley et al.**, Ribosomal multilocus sequence typing: universal characterization of bacteria from domain to strain, K. A. Jolley, C. M. Bliss, J .S. Bennett, H. B. Bratcher, C. Brehony, F. M. Colles, H. Wimalarathna, O. B. Harrison, S. K. Sheppard, A. J. Cody, M .C. Maiden, Microbiology, 2012, [paper link](http://www.ncbi.nlm.nih.gov/pubmed/22282518)
- **Orlek et al.**: Ordering the mob: Insights into replicon and MOB typing schemes from analysis of a curated dataset of publicly available plasmids, A. Orlek, H. Phan, A. E. Sheppard, M. Doumith, M. Ellington, T. Peto, D. Crook, A. S. Walker, N. Woodford, M. F. Anjum, N. Stoesser, Plasmid, 2017, [paper link](https://www.ncbi.nlm.nih.gov/pubmed/28286183)
- **Yutin et al.**: Distribution of ribosomal protein genes across bacterial genome partitions, N. Yutin, P. Puigbò, E. V. Koonin, Y. I. Wolf, PLoS One, 2012, [paper link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3353972/)
- **Villa et al.**: Replicon sequence typing of IncF plasmids carrying virulence and resistance determinants, L. Villa, A. García-Fernández, D. Fortini, A. Carattoli, Journal of Antimicrobial Chemotherapy, 2010, [paper link](https://www.ncbi.nlm.nih.gov/pubmed/20935300)
- **CGE core module**: [repository link](https://bitbucket.org/genomicepidemiology/cge_core_module)
- **fuzzywuzzy**:  [repository link](https://github.com/seatgeek/fuzzywuzzy)
- **MOB-suite**:  “MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies.” Robertson, James, and John H E Nash.  Microbial genomics vol. 4,8 (2018) [paper link](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000206), [repository link](https://github.com/phac-nml/mob-suite)
- **taxize**: "taxize: taxonomic search and retrieval in R." Chamberlain SA and Szöcs E. F1000Res. 2013;2:191. [paper link](https://f1000research.com/articles/2-191/v2), [repository link](https://github.com/ropensci/taxize)
## Notes

This data processing pipeline makes use of the [PubMLST website](https://pubmlst.org/) developed by Keith Jolley ([Jolley & Maiden 2010, BMC Bioinformatics, 11:595](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-595)) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust.
