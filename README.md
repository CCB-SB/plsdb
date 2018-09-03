# Retrieving and processing plasmids from NCBI

## Requirements

### Python modules with Conda

*Requires Python 3*

#### Miniconda
```bash
# get miniconda (for linux)
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# install
bash Miniconda3-latest-Linux-x86_64.sh
# set path to binaries
export PATH=$HOME/miniconda3/bin:$PATH
```

#### Environment and modules
```bash
# create env
conda create --name plsdb python=3
# install packages
conda install --name plsdb -c anaconda -c bioconda -c conda-forge --file requirements.txt
# activate env
source activate plsdb
```

### R packages

Only the R-packages imported in `create_plots.R` are quired. Use `install.packages()` to install missing packages.

### Other software

The binaries of [edirect/eutils](https://www.ncbi.nlm.nih.gov/books/NBK179288/), [Mash](https://github.com/marbl/Mash), and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) are downloaded by the pipeline.

pMLST data is downloaded from [PubMLST](https://pubmlst.org/plasmid/).

### GoogleMaps API key

To map locations of associated BioSamples the GoogleMaps API is used.
Thus, a Google API key is needed.
The key should be stored in a local file specified in the pipeline config (`pipeline.json`, see the entry for `misc/gmaps_api_keys`).
Also, a file with some of the already retrieved locations is included (`locs.pck`, Python 3 `pickle` object) and
will be updated with newly retrieved locations if you run the pipeline.

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

*Note: All relevant files will be created for each source separately.*

- NCBI nucleotide database sources:
    - INSDC (DDBJ, EMBL/ENA, GenBank): [International Nucleotide Sequence Database Collaboration](https://www.ncbi.nlm.nih.gov/genbank/collab/)
    - RefSeq
- Tools:
    - Install [Mash](https://github.com/marbl/Mash)
    - Install [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
    - Install [edirect/eutils](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
    - Get and process [pMLST data from PubMLST DB](https://pubmlst.org/plasmid/)
- Plasmid records:
    - Query for plasmids in the NCBI nucleotide database
    - XX Orlek et al.
- Filtering
    1. XX
    2. XX
        - Retrieve linked assembly UIDs
        - Get assembly level, status, and if it is the latest version
            - [NCBI assembly help page](https://www.ncbi.nlm.nih.gov/assembly/help/)
        - Filtering:
            - No assembly: Completeness status of the nuccore record has to be "complete"
            - Has assembly: assembly status of the latest version has to be "Complete genome"
- BioSamples:
    - Query for BioSamples associated with retrieved plasmids
    - For each BioSample retrieve additional information
        - Process location information
            - Uses GoogleMaps API
            - Uses coordinates if available, otherwise location
        - Note that it can happen that some bio-sample UIDs have no hits which will be printed during the rule execution
- Taxa:
    - Query for taxa associated with retrieved plasmids
    - Process retrieved data to extract queried taxon (ID, name, rank), complete lineage, and taxa/IDs for relevant ranks (from species to superkingdom)
- Plasmid nucleotide sequences:
    - Download FASTA for each plasmid record
    - Create BLAST database file from plasmid FASTA
    - Create sketches from plasmid FASTA using Mash
    - Embedding:
        - Compute pairwise distances between plasmids using Mash
        - Compute embedding using UMAP
            - Requires ca. 42Gb for ca. 37.6K sequences
    - Annotate using ABRicate:
        - BLASTn search in DBs provided by ABRicate
        - Hits are processed and filtered, and collected in one file
    - Annotate using pMLST:
        - Use `mlst` to run BLAST search on downloaded pMLST profiles
        - Process the results
- Create info table:
    - Record information
        - Sequence length and GC content
        - Taxonomic information
        - Other
    - BioSample information
        - Location (as given in DB) and coordinates retrieved with GoogleMaps
        - Isolation source
    - Embedding coordinates
    - PlasmidFinder hits
    - pMLST hits

#### Remarks

**pMLST**

As pMLST is not yet supported by `mlst` the data needs to be dowloaded and pre-processed before it can be used by the tool.
However, some things need to be considered:
- Scheme name/ID:
    - `http://rest.pubmlst.org/db/pubmlst_plasmid_seqdef/classification_schemes` gives an empty list (2018.08.03), i.e. available IDs are unknown
    - `pipeline.json` contains a mapping for scheme names and IDs and it should be checked before running the pipeline
- No profiles:
    - A scheme may have no profiles (e.g. IncF) but `pmlst` requires a non-empty file
    - Thus, a dummy profile needs to be created and the hits to this profile need to be processed accordingly (i.e. rm the dummy ST)
- Problematic profile file formatting
    - `mlst` requires to have an ST column with numeric values and one column per locus
    - E.g. IncA/C cgMLST has "cgST" instead of "ST" and STs in the format "number.number"
    - In such cases the column is renamed and STs are mapped to 1..N (the original values are saved in a diff. column)
    - Here, the results need to be processed to map the ST back to the original ST value

# References

- **Mash**: "Mash: fast genome and metagenome distance estimation using MinHash", B. D. Ondov, T. J. Treangen, P. Melste
d, A. B. Mallonee, N. H. Bergman, S. Koren and A. M. Phillippy, Genome Biology, 2016, [paper link](https://genomebiology
.biomedcentral.com/articles/10.1186/s13059-016-0997-x), [repository link](https://github.com/marbl/Mash)
- **UMAP**: "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction", L. McInnes and J. Healy, arXi
v, 2018,
[paper link](https://arxiv.org/abs/1802.03426), [repository link](https://github.com/lmcinnes/umap)
- **BLAST**: "Basic local alignment search tool." , Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J., J.
 Mol. Biol. 215:403-410, [BLAST paper link](https://www.ncbi.nlm.nih.gov/pubmed/2231712?dopt=Citation), [BLAST+ paper link](https://www.ncbi.nlm.nih.gov/pubmed/20003500), [tool link](https://bl
ast.ncbi.nlm.nih.gov/Blast.cgi)
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
