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

### Other software

The binaries of [edirect/eutils](https://www.ncbi.nlm.nih.gov/books/NBK179288/), [Mash](https://github.com/marbl/Mash), and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) are downloaded by the pipeline.

### GoogleMaps API key

To map locations of associated BioSamples the GoogleMaps API is used.
Thus, a Google API key is needed.
The key should be stored in a local file specified in the pipeline config (`pipeline.json`, see the entry for `misc/gmaps_api_keys`).
Also, a file with some of the already retrieved locations is included (`locs.pck`, Python `pickle` object) and
will be updated with newly retrieved locations if you run the pipeline.

## Pipeline

To print all rules to be executes run:

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
    - EMBL
    - INSDC (DDBJ, ENA, GenBank): [International Nucleotide Sequence Database Collaboration](https://www.ncbi.nlm.nih.gov/genbank/collab/)
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

#### Notes to plasmid query in `nuccore`

- `-molecule genomic` in `efilter` can remove some plasmid sequences
    - Example: NZ_CP013186
    - Thus, not used in the query

#### Plot embedding in R

```R
# load ggplot2 for plotting
require(ggplot2)

# data tag
tag <- '2018_07_25' # replace by your tag

# read in
d_g <- read.csv(file=sprintf('data/master/%s__genbank.umap', tag), header=TRUE, sep='\t')
d_r <- read.csv(file=sprintf('data/master/%s__refseq.umap', tag), header=TRUE, sep='\t')

# plot
pdf('embedding.pdf')
ggplot(data=d_g, aes(x=D1, y=D2)) + geom_point(size=1, alpha=0.5, colour='white', fill='#3399FF', shape=21) + ggtitle('GenBank') + theme_bw()
ggplot(data=d_r, aes(x=D1, y=D2)) + geom_point(size=1, alpha=0.5, colour='white', fill='#3399FF', shape=21) + ggtitle('RefSeq') + theme_bw()
dev.off()
```

# References

- **Mash**: "Mash: fast genome and metagenome distance estimation using MinHash", B. D. Ondov, T. J. Treangen, P. Melsted, A. B. Mallonee, N. H. Bergman, S. Koren and A. M. Phillippy, Genome Biology, 2016, [paper link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x), [repository link](https://github.com/marbl/Mash)
- **UMAP**: "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction", L. McInnes and J. Healy, arXiv, 2018,
[
