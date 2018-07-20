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
        - About chormids: https://www.sciencedirect.com/science/article/pii/S0966842X09002698
        - Plasmid: "plasmid" in FASTA header or unknown type and not the first record
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
tag='2018_07_18'
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

### Embedding
Visuialize in R:

```R
# load ggplot2 for plotting
require(ggplot2)

# path to file containing embedding
tag <- '2018_07_18' # replace by your tag
f <- sprintf('data/master/%s__master_genomic.fna.umap', tag)

# read in
d <- read.csv(file=f, header=TRUE, sep='\t')

# plot
ggplot(data=d, aes(x=D1, y=D2)) + geom_point(size=1, alpha=0.5, colour='white', fill='#3399FF', shape=21) + theme_bw()
```

# EUTILS

Installation:

```bash
perl -MNet::FTP -e \
    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
     $ftp->login; $ftp->binary;
     $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
gunzip -c edirect.tar.gz | tar xf -
rm edirect.tar.gz
cd XX && wget https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz && gunzip -f xtract.Linux.gz && chmod +x xtract.Linux
```

Query:

```bash
time (tools/edirect/esearch -db nuccore -query "plasmid AND "Bacteria"[porgn:__txid2]" | tools/edirect/efetch -format docsum | tools/edirect/xtract -pattern DocumentSummary -element Id Caption Title CreateDate TaxId Biomol MolType Topology SourceDb ProjectId Genome SubName AssemblyAcc Tech Completeness BioSample Organism > test.tmp)
```

https://ncbiinsights.ncbi.nlm.nih.gov/2013/02/19/how-to-download-bacterial-genomes-using-the-entrez-api/

ID -> FASTA: https://www.biostars.org/p/145330/

```bash
#!/bin/bash
split -l 500 -a 5 $1 input.
for f in input.*
do
IDs=$(cat $f | tr "\n" ",")
epost -db nuccore -id $IDs | efetch -format fasta > $f.output
done
cat *.output > plasmids.fna
#Everything worked out? rm input.* output.*
```

https://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=download&orgn=%22%22[orgn]&status=50|40|30|20&report=plasmids&king=All&group=All&subgroup=All&format=undefined

https://www.biostars.org/p/214355/

https://github.com/haruosuz/plasmids/blob/master/scripts/my_data.sh

# Some notes on data retrieval

## BioProject

## BioSample

`Bio.Entrez` can give wrong results, e.g.

```python
from Bio import Entrez
Entrez.email = "Your.Name.Here@example.org"
handle = Entrez.efetch(db="biosample", id='SAMEA2272644', retmode="xml")
handle.readlines()
# will have: accession="SAMN02272644" which is a different sample
```

# Getting plasmids from NCBI

Using `eutils` for queries ([reference](https://www.ncbi.nlm.nih.gov/books/NBK179288/)).

Why:

- Using assembly reports
    - Many plasmids are submitted as already fully assembled sequences, thus they are not included in the reports
- Using RefSeq plasmid [report](ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt) and [FASTAs](ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/)
    - Not all sequences from the FASTAs are in the table and vice versa
    - Only available for RefSeq, many GenBank plasmids without a RefSeq ID are not in the report
- `Entrez` module of the `Biopython` library
    - Queries do not always work correctly, e.g. getting results for a wrong ID

## Plasmid queries

Reference: https://www.ncbi.nlm.nih.gov/books/NBK179288/
- Parameters:
    - e.g. for efilter: "Use efilter to restrict search or link results by indexed terms"
    - e.g. for einfo: "The einfo function returns information on Entrez indexed fields:"

```bash
./esearch -db nuccore -query "plasmid" | ./efilter -location plasmid -molecule genomic -organism bacteria -source refseq | ./efetch -format docsum |./xtract -pattern DocumentSummary -ID "(NA)" -CAP "(NA)" -TITLE "(NA)" -DATE "(NA)" -TOPO "(NA)" -COMPL "(NA)" -PID "(NA)" -SID "(NA)" -AID "(NA)" -TID "(NA)" -GEN "(NA)" -ID Id -CAP Caption -TITLE Title -DATE CreateDate -TOPO Topology -COMPL Completeness -PID ProjectId -SID BioSample -AID AssemblyAcc -TID TaxId -GEN Genome -tab "\t" -element "&ID" "&CAP" "&TITLE" "&DATE" "&TOPO" "&COMPL" "&PID" "&SID" "&AID" "&TID" "&GEN"
```

## Plasmid FASTAs

Reference: https://ncbiinsights.ncbi.nlm.nih.gov/2013/02/19/how-to-download-bacterial-genomes-using-the-entrez-api/

```bash
epost -db nuccore -id $IDs | efetch -format fasta > $f.output
```

## BioSample

```bash
./epost -db biosample -id SAMN09388963 -format acc | ./efetch -format xml | ./xtract -pattern BioSample  -LOC "(NA)" -COORD "(NA)" -SOURCE "(NA)" -ACC @accession -block Attribute -if Attribute@attribute_name -equals 'geo_loc_name' -LOC Attribute -block Attribute -if Attribute@attribute_name -equals 'lat_lon' -COORD Attribute -block Attribute -if Attribute@attribute_name -equals 'isolation_source' -SOURCE Attribute -tab "\t" -element "&ACC" "&LOC" "&COORD" "&SOURCE"
```

## Taxonomy

he "*/Child" construct will skip past the outer start tag:

  efetch -db taxonomy -id 9606,7227,10090 -format xml |
  xtract -pattern Taxon -block "*/Taxon" \
    -tab "\n" -element TaxId,ScientificName
to visit the next level of nested objects individually:

  131567    cellular organisms
  2759      Eukaryota
  33154     Opisthokonta
  ...
