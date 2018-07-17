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
# check
comm -3 <(pip freeze | sort) <(sort requirements.txt)
```

### Other software

The binaries of [Mash](https://github.com/marbl/Mash) and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) are downloaded by the pipeline.

### GoogleMaps API queries

To map locations of associated BioSamples the GoogleMaps API is used.
Thus, a Google API key is needed.
The keys should be stored in a local file specified in the pipeline config
(`pipeline.json`, see the keys `misc`, `gmaps_api_keys`).
Also, as there are many location one run is not sufficient
due to the max. limit of queries per day set by Google.
Thus, a file with some of the already retrieved locations is included (`locs.pck`).
Otherwise, the script for creating the master table should be run multiple times on different days.

## Pipeline

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
    - Mash distances from the Mash sketch
    - Embedding using UMAP from Mash distances
    - File with all genomic features
    - Table with meta data
        - Assembly information
        - Taxonomic information from assembly taxon
        - Sequence information (length, GC)
        - BioSample information

## Output checks

### Downloaded genomic sequences
To check if all required FASTAs were successfully downloaded and processed.

**Count downloads**
```bash
# replace by your tag, i.e. date (year_month_day) when the data was downloaded
tag='2018_06_28'
f="data/genomes/${tag}/download.log"
# numbre of FASTAs in filtered report files
refseq=$(tail -n +2 data/reports/${tag}__refseq__assembly_report__filt.tsv | wc -l)
genbank=$(tail -n +2 data/reports/${tag}__genbank__assembly_report__filt.tsv | wc -l)

# stats
# expected number of executed download CMDs
echo "expected: $(( ${refseq} + ${genbank} ))"
# all files to be downloaded
# -> should be = expected
echo "total: $(grep 'CMD: curl' ${f} | wc -l)"
# failed downloads
# -> should be 0
echo "failed: $(grep 'FAILED DOWNLOAD' ${f} | wc -l)"
# server file does not exist
# -> should be 0
echo "no file: $(grep 'NO SUCH FILE' ${f} | wc -l)"
# local files already existed so not downloaded
# -> should be 0
echo "existed: $(grep 'ALREADY EXISTS' ${f} | wc -l)"
# could not process FASTA to extract plasmids
# -> should be 0
echo "error: $(grep 'ERROR' ${f} | wc -l)"
```

**Check if plasmids were extracted from all files**
```bash
# FASTAs from wget CMDs
grep 'CMD: curl' ${f} | grep -o '\-\-output .* ftp://' | cut -d' ' -f2 | sort | uniq > check_fastas_all.tmp
# FASTAs without plasmids
grep -o 'NO PLASMIDS.*$' ${f} | cut -d' ' -f4 | sort | uniq > check_fastas_no.tmp
# FASTAs with plasmids
grep 'INFO - PLASMID' ${f} | grep -o 'data/genomes/.*$' | cut -d' ' -f3 | sort | uniq > check_fastas_pls.tmp
# FASTAs which should have plasmids
comm -3 check_fastas_all.tmp check_fastas_no.tmp > check_fastas_should.tmp
# FASTAs with plasmids but not saved
comm -3 check_fastas_should.tmp check_fastas_pls.tmp > check_fastas_miss.tmp
# -> should be 0
echo "$(wc -l check_fastas_miss.tmp)"
# count files
# -> should be (all - no plasmids)
echo "files: $(find data/genomes/${tag}/ -type f -name '*_genomic.fna.gz' | wc -l)"
```

### Downloaded genomic features
Files: `*_feature_table.txt.gz`

For some assemblies this file may be missing thus trying to download it will fail.
Check the number of successful and failed download attempts.

```bash
tag='2018_06_28'
f="data/features/${tag}/add_download.log"

# all files to be downloaded
# -> should be = total number of FASTAs in filtered reports
echo "total: $(grep 'CMD: curl' ${f} | wc -l)"
# failed downloads
# -> should be 0
echo "failed: $(grep 'FAILED DOWNLOAD' ${f} | wc -l)"
# server file does not exist
# -> may be > 0
echo "no file: $(grep 'NO SUCH FILE' ${f} | wc -l)"
# local files already existed so not downloaded
# -> should be 0
echo "existed: $(grep 'ALREADY EXISTS' ${f} | wc -l)"
# could not process FASTA to extract plasmids
# -> should be 0
echo "error: $(grep 'ERROR' ${f} | wc -l)"

# count downloaded files
# -> should be (all - no file)
echo "files: $(find data/features/${tag}/ -type f -name '*_feature_table.txt.gz' | wc -l)"
```

### Embedding
Visuialize in R:

```R
# load ggplot2 for plotting
require(ggplot2)

# path to file containing embedding
tag <- '2018_06_28' # replace by your tag
f <- sprintf('data/master/%s__master_genomic.fna.umap', tag)

# read in
d <- read.csv(file=f, header=TRUE, sep='\t')

# plot
ggplot(data=d, aes(x=D1, y=D2)) + geom_point(size=1, alpha=0.5, colour='white', fill='#3399FF', shape=21) + theme_bw()
```
