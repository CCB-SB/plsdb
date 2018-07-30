# Retrieving and processing plasmids from NCBI

## Requirements

### Python modules
Create a python environment and install required modules:

```bash
# create and activate venv
virtualenv -p /usr/bin/python3 venv
source venv/bin/activate
# requirements
pip install -r requirements.txt
```

### R packages

Only the R-packages imported in `create_plots.R` are quired. Use `install.packages()` to install missing packages.

### Other software

The binaries of [edirect/eutils](https://www.ncbi.nlm.nih.gov/books/NBK179288/), [Mash](https://github.com/marbl/Mash), and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) are downloaded by the pipeline.

### GoogleMaps API key

To map locations of associated BioSamples the GoogleMaps API is used.
Thus, a Google API key is needed.
The key should be stored in a local file specified in the pipeline config (`pipeline.json`, see the entry for `misc/gmaps_api_keys`).
Also, a file with some of the already retrieved locations is included (`locs.pck`, Python `pickle` object) and
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
- Plasmid records:
    - Query for plasmids in the NCBI nucleotide database
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
- Create an info table containing:
    - Record information
        - Sequence length and GC content
        - Taxonomic information
        - Other
    - BioSample information
        - Location (as given in DB) and coordinates retrieved with GoogleMaps
        - Isolation source
    - Embedding coordinates

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
