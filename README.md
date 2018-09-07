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

## Pipeline
To print all rules to be executed run:

```bash
snakemake -s pipeline.snake -np
```

Call the pipeline using
```bash
snakemake -s pipeline.snake
```

### Pipeline steps
- NCBI nucleotide database sources:
    - INSDC (DDBJ, EMBL/ENA, GenBank): [International Nucleotide Sequence Database Collaboration](https://www.ncbi.nlm.nih.gov/genbank/collab/)
    - RefSeq
- Tools:
    - Install Mash
    - Install BLAST+
    - Install edirect/eutils
    - Get and process pMLST data from PubMLST DB
    - Update data for ABRicate
- Plasmid records:
    - Query for plasmids in the NCBI nucleotide database
        - `esearch` query from Orlek *et al.*
- Plasmid meta data
    - Retrieve linked assemblies and relevant meta data
    - Retrieve linked BioSamples and relevant meta data
    - Retrieve taxonomic information
        - Process: extract queried taxon (ID, name, rank), complete lineage, and taxa/IDs for relevant ranks (from species to superkingdom)
    - Add all new meta data to the table
    - Process location information of the BioSamples and add it to the table
        - Use coordinates if available, otherwise location
        - Use OpenCageData API
- Filtering (1): Filter by
    - Record description (regular expression from Orlek *et al.*)
    - (Assembly) completeness
        - If no assembly: Completeness status of the nuccore record has to be `complete`
        - Has assembly: assembly status of the latest version has to be `Complete genome`
            - [NCBI assembly help page](https://www.ncbi.nlm.nih.gov/assembly/help/)
    - By taxonomy: superkingdom taxon ID should be `2` (i.e. Bacteria)
- Filtering (2): Remove identical plasmids
    - Download nucl. sequences of plasmid records
    - Compute the sketches using Mash
    - Get pairs of plasmids with distance of 0 using Mash
    - Group plasmids with identical sequences
    - Among these groups select one record
        - Prefer RefSeq records and those with more information
- Filtering (3): Find and remove putative chromosomal sequences
    - Create the BlastDB of rMLST allele sequences
        - The FASTAs are **NOT** downloaded by the pipeline (see the "Requirements" section above)
    - The plasmid sequences are aligned against the rMLST allele sequences
    - Remove plasmids having more than 5 unique rMLST loci
- Plasmid nucleotide sequences:
    - Create a new FASTA with nucl. sequences of remained plasmids
    - Annotate using ABRicate:
        - BLASTn search in DBs provided by ABRicate
        - Hits are processed and filtered, and collected in one file
    - Annotate using pMLST:
        - Use `mlst` to run BLAST search on downloaded pMLST profiles
        - Process the results
    - Create BLAST database file from plasmid FASTA
    - Create sketches from plasmid FASTA using Mash
- Embedding:
    - Compute pairwise distances between plasmids using Mash
    - Compute embedding using UMAP
        - Requires ca. XXGb for ca. XXK sequences
- Create info table:
    - Record information
    - Embedding coordinates
    - PlasmidFinder hits
    - pMLST hits

#### Notes

##### rMLST

[XX](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3353972/): "Distribution of ribosomal protein genes across bacterial genome partitions"

    In 68 of the 995 analyzed bacterial genomes, r-protein genes are
    distributed across two or more genome partitions. In some cases,
    paralogous proteins are encoded in different chromosomes or plasmids.
    [...]
    In other cases, r-protein genes are present in a single copy
    but are spread across genome partitions.

##### pMLST
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
- **PlasmidFinder**: "In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing.", A. Carattoli, E. Zankari, A. Garcia-Fernandez, M. Voldby Larsen, O. Lund, L. Villa, F. Møller Aarestrup, and H. Hasman, Antimicrob. Agents Chemother., 2014, [paper link](http://aac.asm.org/content/58/7/3895.long)
- **pMLST in PubMLST**: [web-site](https://pubmlst.org/plasmid/)
- **mlst**: Tool implemented by Thorsten Seemann, [repository link](https://github.com/tseemann/mlst)
- **checkM**: CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes, D. H. Parks, M. Imelfort, C. T. Skennerton, P. Hugenholtz and G. W. Tyson, Genome Res., 2015, [paper link](https://genome.cshlp.org/content/25/7/1043.short), [repository link](https://github.com/Ecogenomics/CheckM)
- **OpenCageData**: An API to convert coordinates to and from places, [web-site](https://opencagedata.com/)
- **rMLST**: [rMLST at PubMLST](https://pubmlst.org/rmlst/), [paper link](http://www.ncbi.nlm.nih.gov/pubmed/22282518)
- **Orlek et al.**: Ordering the mob: Insights into replicon and MOB typing schemes from analysis of a curated dataset of publicly available plasmids, A. Orlek, H. Phan, A. E. Sheppard, M. Doumith, M. Ellington, T. Peto, D. Crook, A. S. Walker, N. Woodford, M. F. Anjum, N. Stoesser, Plasmid, 2017, [paper link](https://www.ncbi.nlm.nih.gov/pubmed/28286183)
