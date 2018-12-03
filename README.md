# Retrieving and processing plasmids from NCBI

## Requirements

### Python

#### Miniconda
```bash
cd ~
# get miniconda (for linux)
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# install
bash Miniconda3-latest-Linux-x86_64.sh
# set path to binaries
export PATH=$HOME/miniconda3/bin:$PATH
```

#### Main conda environment
Create the main environment and install needed packages:
```bash
# create env
conda create --name plsdb python=3
# install packages
conda install --name plsdb -c anaconda -c bioconda -c conda-forge --file requirements.txt
# activate env
source activate plsdb
```

**IMPORTANT**: Currently, ABRicate (version `0.8.10`) does not download the PlasmidFinder sequences from the BitBucket repository ([see this issue](https://github.com/tseemann/abricate/issues/66)).
Thus, the code in the function `get_plasmidfinder` in `abricate-get_db` has to be changed to download the most recent version of the sequence files.

```bash
rsync -av patch_abricate-get_db $HOME/miniconda3/envs/plsdb/bin/abricate-get_db
```

### R packages
Only the R-packages imported in `create_plots.R` are quired.
Use `install.packages()` to install missing packages.

### Other tools
Other tools installed by the pipeline:
- Binaries of [edirect/eutils](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- Binaries of [Mash](https://github.com/marbl/Mash)
- Binaries of [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+)

### Datasets
The pipeline will automatically update the following databases/datasets:
- pMLST data is downloaded from [PubMLST](https://pubmlst.org/plasmid/)
- ABRicate data is updated using the built-in function

**IMPORTANT**: The rMLST sequences from PubMLST are **NOT** downloaded by the pipeline as the access to the sequences requires a login.
The pipeline expects a single FASTA file with all sequences (its path should be set in the config file `pipeline.json`, see `rmlst/fas`).

### API key for location queries
To map location names to coordinates the [OpenCageData](https://opencagedata.com/) API is used which requires an API key (you can register for a free trial account).
The key should be stored in a local file specified in the pipeline config (`pipeline.json`, see `data/api_keys`).
Also, a file with some of the already retrieved locations is included (`locs.tsv`) and will be updated with newly retrieved locations if you run the pipeline.

--------------------------------------------------

## Pipeline
To print all rules to be executed run:

```bash
snakemake -s pipeline.snake -np
```

Call the pipeline using
```bash
snakemake -s pipeline.snake
```

**IMPORTANT**: It is better to run the pipeline step by step to perform manual checking of the created files and the logger output.

### Steps
- Used NCBI nucleotide database sources:
    - INSDC (DDBJ, EMBL/ENA, GenBank): [International Nucleotide Sequence Database Collaboration](https://www.ncbi.nlm.nih.gov/genbank/collab/)
    - RefSeq
- Tools/data:
    - Install Mash
    - Install BLAST+
    - Install edirect/eutils
    - Get and process pMLST data from PubMLST DB
    - Update data for ABRicate
- Plasmid records:
    - Query for plasmids in the NCBI nucleotide database
        - `esearch` query from **Orlek et al.**
- Plasmid meta data
    - Retrieve linked assemblies and relevant meta data
    - Retrieve linked BioSamples and relevant meta data
    - Retrieve taxonomic information
        - Process: extract queried taxon (ID, name, rank), complete lineage, and taxa/IDs for relevant ranks (from species to superkingdom)
    - Add all new meta data to the table
    - Process location information of the BioSamples and add it to the table
        - Use coordinates if available, otherwise location
        - Use OpenCageData API
- Filtering (1): To remove incomplete or nonbacterial records these are filtered by
    - Record description (regular expression from **Orlek et al.**)
    - (Assembly) completeness
        - If no assembly: Completeness status of the nuccore record has to be `complete`
        - Has assembly: assembly status of the latest version has to be `Complete genome`
            - [NCBI assembly help page](https://www.ncbi.nlm.nih.gov/assembly/help/)
    - By taxonomy: superkingdom taxon ID should be `2` (i.e. Bacteria)
- Filtering (2): To remove identical records
    - Download nucl. sequences of plasmid records
    - Compute the sketches using Mash
    - Get pairs of plasmids with distance of 0 using Mash
    - Group plasmids with identical sequences
    - Among these groups select one record
        - Prefer RefSeq records and those with more information
- Filtering (3): To remove putative chromosomal sequences
    - Create the BlastDB of rMLST allele sequences
        - The FASTAs are **NOT** downloaded by the pipeline (see the "Requirements" section above)
    - The plasmid sequences are aligned against the rMLST allele sequences
    - Records having more than 5 unique rMLST loci are searched in NCBI chromosomal sequences using BLASTn (remote access)
    - Records with hits are removed
- Plasmid nucleotide sequences:
    - Create a new FASTA with nucl. sequences of remained plasmids
    - Annotate using ABRicate:
        - BLASTn search in DBs provided by ABRicate
            - Blaster from [CGE core module](https://bitbucket.org/genomicepidemiology/cge_core_module) is used for search and pre-processing
        - Filtering:
            - Identity and coverage cutoffs
            - Overlapping matches are removed
        - All hits are ollected into one file
    - Annotate using pMLST:
        - For each found replicon use the associated pMLST scheme (if available)
        - Use `mlst` to perform the pMLST analysis
        - Process the results
            - Set IncF ST according to the FAB formula (**Villa et al.**)
    - Create BLAST database file from plasmid FASTA
    - Create sketches from plasmid FASTA using Mash
- List of similar plasmids:
    - Use Mash do compute pairwise distances (use a distance cutoff)
    - Create a list of unique pairs
- Embedding:
    - Compute pairwise distances between plasmids using Mash
    - Compute embedding using UMAP
- Create info table:
    - Record information
    - Embedding coordinates
    - PlasmidFinder hits
    - pMLST hits

#### Notes

##### Finding putative chromosomal sequences

The candidates for putative chromosomal sequences are determined by searching for the rps genes - ribosome protein subunits which are used in the rMLST scheme (containing 53 rps genes) introduced by **Jolley et al.**:

    "The rps loci are ideal targets for a universal characterization scheme as they are:
      (i) present in all bacteria;
     (ii) distributed around the chromosome; and
    (iii) encode proteins which are under stabilizing selection for functional conservation."

However, some of the rps genes can also be found on plasmids as described by **Yutin et al.**:

    In 68 of the 995 analyzed bacterial genomes, r-protein genes are
    distributed across two or more genome partitions. In some cases,
    paralogous proteins are encoded in different chromosomes or plasmids.

Thus, the presence of (some) rps genes alone cannot be always used as an indicator for chromosomal sequences.
Therefore, the records containing more than 5 rps genes are searched in the NCBI sequences using BLAST.

##### Performing pMLST with the tool "mlst"
As pMLST is not yet supported by `mlst` the data needs to be dowloaded and pre-processed before it can be used by the tool.
However, some things need to be considered:
- No profiles:
    - A scheme may have no profiles (e.g. IncF) but `pmlst` requires a non-empty file
    - Thus, a dummy profile needs to be created and the hits to this profile need to be processed accordingly (i.e. by removing the dummy ST)
- Problematic profile file formatting:
    - `mlst` requires an ST column with numeric values and one column per locus
    - E.g. IncA/C cgMLST has "cgST" instead of "ST" and STs in the format "number.number"
    - In such cases the column is renamed and STs are mapped to 1..N (the original values are saved in a different file)
    - Here, the results need to be processed to map the ST back to the original ST value

Also make sure to provide a mapping for each downloaded pMLST scheme to a Python regular expression to be used to match plasmidFinder hits and the pMLST scheme name (see `pmlst/map` in `pipeline.json`).

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
- **checkM**: CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes, D. H. Parks, M. Imelfort, C. T. Skennerton, P. Hugenholtz and G. W. Tyson, Genome Res., 2015, [paper link](https://genome.cshlp.org/content/25/7/1043.short), [repository link](https://github.com/Ecogenomics/CheckM)
- **OpenCageData**: An API to convert coordinates to and from places, [web-site](https://opencagedata.com/)
- **rMLST**: [rMLST at PubMLST](https://pubmlst.org/rmlst/)
- **Jolley et al.**, Ribosomal multilocus sequence typing: universal characterization of bacteria from domain to strain, K. A. Jolley, C. M. Bliss, J .S. Bennett, H. B. Bratcher, C. Brehony, F. M. Colles, H. Wimalarathna, O. B. Harrison, S. K. Sheppard, A. J. Cody, M .C. Maiden, Microbiology, 2012, [paper link](http://www.ncbi.nlm.nih.gov/pubmed/22282518)
- **Orlek et al.**: Ordering the mob: Insights into replicon and MOB typing schemes from analysis of a curated dataset of publicly available plasmids, A. Orlek, H. Phan, A. E. Sheppard, M. Doumith, M. Ellington, T. Peto, D. Crook, A. S. Walker, N. Woodford, M. F. Anjum, N. Stoesser, Plasmid, 2017, [paper link](https://www.ncbi.nlm.nih.gov/pubmed/28286183)
- **Yutin et al.**: Distribution of ribosomal protein genes across bacterial genome partitions, N. Yutin, P. Puigbò, E. V. Koonin, Y. I. Wolf, PLoS One, 2012, [paper link]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3353972/)
- **Villa et al.**: Replicon sequence typing of IncF plasmids carrying virulence and resistance determinants, L. Villa, A. García-Fernández, D. Fortini, A. Carattoli, Journal of Antimicrobial Chemotherapy, 2010, [paper link](https://www.ncbi.nlm.nih.gov/pubmed/20935300)
- **CGE core module**: [repository link](https://bitbucket.org/genomicepidemiology/cge_core_module)

The data processing pipeline makes use of the [PubMLST website](https://pubmlst.org/) developed by Keith Jolley ([Jolley & Maiden 2010, BMC Bioinformatics, 11:595](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-595)) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust.
