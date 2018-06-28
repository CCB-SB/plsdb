# Retrieving and processing plasmids from NCBI

## Pipeline summary

Create a python environment and install required modules:

```bash
# create and activate venv
virtualenv -p /usr/bin/python3 venv
source venv/bin/activate
# other requirements
pip install -r requirements.txt
# check
comm -3 <(pip freeze | sort) <(sort requirements.txt)
```

Call the pipeline using `snakemake -s pipeline.snake`.

### Pipeline steps

- Preliminary steps:
    - Install [Mash](https://github.com/marbl/Mash)
    - Install [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- Download assembly reports (RefSeq, GenBank)
    - [Bacteria assembly summary, GenBank](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt)
    - [Bacteria assembly summary, RefSeq](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt)
- Filter assembly reports
    - Status `latest`
    - No scaffolds/contigs
- Download plasmid FASTAs: `*_genomic.fna.gz`
    - Download FASTA files with nucl. sequences
    - Filter out non-plasmid sequences
- Create table of downloaded plasmid sequences
- Download and filter other assembly data
    - Genomic features: `*_feature_table.txt.gz`
- Create "master" files
    - FASTA with all nucl. sequences
    - BLAST DB from the FASTA file
    - Mash sketch from the FASTA file
    - Mash distances
    - Embedding using UMAP
    - File with all genomic features
    - Table with meta data
        - Assembly information
        - Taxonomic information
        - Sequence information (length, GC)
        - Count genomic features

## Output checks

### Downloaded genomic sequences
To check of all required FASTAs were successfully downloaded and processed.

```bash
tag='2018_06_25'
f="data/genomes/${tag}/download.log"

# FASTAs from wget CMDs
grep 'CMD: wget' ${f} | grep -o '\-O .* ftp://' | cut -d' ' -f2 | sort | uniq > check_fastas_all.tmp
# -> should be the number of FASTAs in filtered RefSeq and Genbank reports
wc -l check_fastas_all.tmp

# Failed downloads
# -> should be 0
grep 'FAILED DOWNLOAD' ${f} | wc -l

# FASTAs without plasmids
grep -o 'NO PLASMIDS.*$' ${f} | cut -d' ' -f4 | sort | uniq > check_fastas_no.tmp
# FASTAs with plasmids
grep 'INFO - PLASMID' ${f} | grep -o 'data/genomes/.*$' | cut -d' ' -f3 | sort | uniq > check_fastas_pls.tmp
# FASTAs wich should have plasmids
comm -3 check_fastas_all.tmp check_fastas_no.tmp > check_fastas_should.tmp

# Number of FASTAs with plasmids but not saved
# -> should be 0
comm -3 check_fastas_should.tmp check_fastas_pls.tmp | wc -l
```

### Downloaded genomic features
Files: `*_feature_table.txt.gz`

For some assemblies these file may be missing thus trying to download it will fail.
Check the number of successful and failed download attempts.

```bash
tag='2018_06_25'
f="data/features/${tag}/add_download.log"

# all wget CMDs
grep 'CMD: wget' ${f} | wc -l
# no plasmids in downloaded file
# -> should be 0
grep 'NO PLASMIDS' ${f} | wc -l
# failed wget CMDs (e.g. no annotated features)
grep 'FAILED' ${f} | wc -l
# count downloaded files
# -> should be (all - failed)
find data/features/${tag}/ -type f -name '*_feature_table.txt.gz' | wc -l
```

### Embedding
Plot in R:

```R
# load ggplot2 for plotting
require(ggplot2)
# path to file containing embedding
f <- '<tag>__master_genomic.fna.umap'
# read in
d <- read.csv(file=f, header=FALSE, sep='\t')
# plot
ggplot(data=d, aes(x=D1, y=D2)) + geom_point(size=1, alpha=0.5, colour='white', fill='#3399FF', shape=21) + theme_bw()
```
