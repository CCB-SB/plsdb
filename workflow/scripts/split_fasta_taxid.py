import pandas as pd
from Bio import SeqIO
from utilsmeta.utils import load_table

# Get TAXONOMY_UID for archaea and bacteria
tax = load_table(snakemake.input.tax, table="TAXONOMY")
bact_uids = tax[tax.TAXONOMY_superkingdom_id == 2].TAXONOMY_UID
arc_uids = tax[tax.TAXONOMY_superkingdom_id == 2157].TAXONOMY_UID

# Get NUUCORE_UIDs fro archae and bacteria
nucc = load_table(snakemake.input.nucc, table="NUCCORE")
bact = list(nucc[nucc.TAXONOMY_UID.isin(bact_uids)].NUCCORE_ACC)
arc = list(nucc[nucc.TAXONOMY_UID.isin(arc_uids)].NUCCORE_ACC)

# Filter fasta
records = SeqIO.to_dict(SeqIO.parse(snakemake.input.fasta, "fasta"))
records_bac = {k:v for k,v in records.items() if k in bact}
records_arc = {k:v for k,v in records.items() if k in arc}

with open(snakemake.output.bac, 'w') as out:
    [SeqIO.write(record, out, 'fasta') for k, record in records_bac.items()]
with open(snakemake.output.arc, 'w') as out:
    [SeqIO.write(record, out, 'fasta') for k, record in records_arc.items()]