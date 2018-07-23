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
The keys should be stored in a local file specified in the pipeline config
(`pipeline.json`, see the entry for `misc/gmaps_api_keys`).
Also, a file with some of the already retrieved locations is included (`locs.pck`, Python `pickle` object) and it will be updated with newly retrieved locations if you run the pipeline.

## Pipeline

First, print all rules to be executed: `snakemake -s pipeline.snake -np`.

Call the pipeline using `snakemake -s pipeline.snake`.

### Pipeline steps

*Note: All relevant files will be created for plasmids from GenBank and RefSeq respectively.*

- Tools:
    - Install [Mash](https://github.com/marbl/Mash)
    - Install [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
    - Install [edirect/eutils](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- Plasmid records:
    - Query for plasmids in the NCBI database
    - Filter plasmids by by sequence length, and `Completeness` and `Genome` tags
- BioSamples:
    - Query for BioSamples associated with retrieved plasmids
    - Process locations
        - Uses GoogleMaps API
        - Location coordinates are used for mapping if available, otherwise location name is used
- Taxa:
    - Query for taxa associated with retrieved plasmids
    - Process retrieved data to exract queried taxon (ID, name, rank), complete lineage, and taxa/IDs for relevant ranks (from species to superkingdom)
- Plasmid sequences:
    - Download FASTA for each plasmid record
    - Create BLAST database file from plasmid FASTA
    - Create sketches from plasmid FASTA using Mash
    - Embedding:
        - Compute pairwise distances between plasmids using Mash
        - Compute embedding using UMAP
            - Requires ca. 35Gb for ca. 36K plasmids
    - Create an info table containing:
        - Record information
            - Sequence length and GC content
            - Accession for BioProject, BioSample, Assembly (if available)
            - Taxonomic information
        - BioSample information
            - Location (as given in DB) and coordinates retrieved with GoogleMaps
            - Isolation source
        - Embedding coordinates

#### Output checks

TODO

##### Embedding - plot in R

```R
# load ggplot2 for plotting
require(ggplot2)

# data tag
tag <- '2018_07_23' # replace by your tag

# read in
d_g <- read.csv(file=sprintf('data/master/%s__genbank.umap', tag), header=TRUE, sep='\t')
d_r <- read.csv(file=sprintf('data/master/%s__refseq.umap', tag), header=TRUE, sep='\t')

# plot
pdf('embedding.pdf')
ggplot(data=d_g, aes(x=D1, y=D2)) + geom_point(size=1, alpha=0.5, colour='white', fill='#3399FF', shape=21) + ggtitle('GenBank') + theme_bw()
ggplot(data=d_r, aes(x=D1, y=D2)) + geom_point(size=1, alpha=0.5, colour='white', fill='#3399FF', shape=21) + ggtitle('RefSeq') + theme_bw()
dev.off()
```
