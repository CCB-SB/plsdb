from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import gc_fraction
import pandas as pd
import numpy as np
from os.path import join
import os
from utilsmeta.utils import load_table
 
# Load Fasta and nuccore
fasta = SeqIO.to_dict(SeqIO.parse(snakemake.input.fasta, "fasta"))
records = SeqIO.to_dict(SeqIO.parse(snakemake.input.genbank, "genbank"))
nucc = load_table(snakemake.input.nucc, table="NUCCORE")
nucc.set_index("NUCCORE_ACC", drop=False, inplace=True)

#
## AMR results
#

df = pd.read_table(snakemake.input.amr, low_memory=False)
df['NUCCORE_ACC'] = df.loc[:, 'input_sequence_id']
df['strand'] = np.where(df['strand_orientation'] == '-', -1, +1)

# Detect equals gene symbols. Example: OXA-48 (rgi) and blaOXA-48 (ncbi)
ncbi_df = df[df['analysis_software_name'] == 'amrfinderplus']
bla_ncbi_genes = set(i.lower().replace('bla', '') for i in set(ncbi_df['gene_symbol']) if i.startswith('bla'))
df['gene_symbol2'] = np.where((df['analysis_software_name'] == 'rgi') & (df['gene_symbol'].str.lower().isin(bla_ncbi_genes)),
         'bla'+df['gene_symbol'].str.lower(), df['gene_symbol'].str.lower())

cols = ["NUCCORE_ACC", 
    "analysis_software_name",
    "reference_database_version",
    "gene_symbol", "gene_symbol2",
    "gene_name",
    "drug_class",
    "antimicrobial_agent",
    "input_gene_start", 
    "input_gene_stop", 
    "input_gene_length",
    "strand_orientation",
    "sequence_identity",
    "coverage_percentage",
    "predicted_phenotype_confidence_level" 
    ]

for index, row in df.iterrows():
    try:
        acc = row['NUCCORE_ACC']
        start = row['input_gene_start']
        stop = row['input_gene_stop']

        if stop < start:
            stop = row['input_gene_start']
            start = row['input_gene_stop']
        
        records[acc].seq = fasta[acc].seq
        
        feature = SeqFeature(
            location = FeatureLocation(start = start, 
                            end = stop,
                            strand = row['strand']), 
            type='AMR', id=row['gene_symbol'],
            qualifiers = {
                "gene": row['gene_symbol'],
                "gene_name": row['gene_name'],
                "tool": row['analysis_software_name'],
                "drug_class": row['drug_class'],
                "identity": row['sequence_identity'],
                "coverage": row['coverage_percentage']}
        )

        ids = {item.qualifiers['gene'][0]: i for i, item in enumerate(records[acc].features) if 'qualifiers' in dir(item) if "gene" in item.qualifiers}
        ids2 = {('bla'+item.qualifiers['gene'][0].lower() if item.qualifiers['gene'][0].lower() in bla_ncbi_genes else item.qualifiers['gene'][0].lower() ): i for i, item in enumerate(records[acc].features) if 'qualifiers' in dir(item) if "gene" in item.qualifiers}
        
        
        if row['gene_symbol'] in ids:
            records[acc].features[ids[row['gene_symbol']]] = feature
        if row['gene_symbol'] in ids2:
            records[acc].features[ids2[row['gene_symbol']]] = feature
        elif row['gene_symbol2'] in ids:
            records[acc].features[ids[row['gene_symbol2']]] = feature
        elif row['gene_symbol2'] in ids2:
            records[acc].features[ids2[row['gene_symbol2']]] = feature
        else:
            if acc == "AJ863570.1" and row['gene_symbol']=="blaOXA-2":
                print("Duplicacted no detected")
                print(row['gene_symbol2'])
                print(ids)
                print(ids2)
                raise Exception

            records[acc].features.append(feature)
                
    except Exception as e:
        print(row)
        raise e


# TABULAR
tab_amr = df.loc[:, cols]
tab_amr.to_csv(snakemake.output.amr_tab, index=False, sep='\t')

#
## ANTISMASH Results
#

# Adjust names
df = pd.read_table(snakemake.input.bgc)
df['NUCCORE_ACC'] = df.loc[:, 'GBK'].replace('\\.region.*\\.gbk', '', regex=True)


for index, row in df.iterrows():
    try:
        acc = row['NUCCORE_ACC']

        if not records[acc].seq.defined:
            records[acc].seq = fasta[acc].seq
        

        feature = SeqFeature(
            location = FeatureLocation(start = row['From'], 
                            end = row['To'],
                            strand = None), 
            type='BGC', id=row['BGC_Type'],
            qualifiers = {
                "gene": row['BGC_Type'],
                "tool": "antismash"}
        )
        records[acc].features.append(feature)

        
    except Exception as e:
        print(row)
        raise e

#
## Save INFO
#
os.makedirs(snakemake.output.DIR)
rows_gc = []
rows_qualifiers = []
protein_seqs = []

qualifiers_name = ['gene', 'locus_tag', 'product', 
                   'protein_id', 'codon_start', 'transl_table', 
                   'GO_process']

with open(snakemake.output.proteins, 'w') as prot_fasta: 
    for acc, record in records.items():
        if not record.seq.defined:
                record.seq = fasta[acc].seq
        
        # Save GC content and Topology
        rows_gc.append((acc,gc_fraction(record.seq), record.annotations['topology']))

        # Save qualifiers info
        if 'features' in dir(record):
            for item in record.features:
                if 'qualifiers' in dir(item):
                    row = [acc]
                    row.extend([item.qualifiers[i][0] if i in item.qualifiers else None for i in qualifiers_name])
                    row.extend([int(item.location.start), int(item.location.end), item.location.strand])
                    rows_qualifiers.append(row)

                    # Save protein sequences
                    try:
                        if 'translation' in item.qualifiers:
                            
                            r = SeqRecord(
                                Seq(item.qualifiers['translation'][0]),
                                id=item.qualifiers['protein_id'][0], description=acc
                                )

                            SeqIO.write(r, prot_fasta, 'fasta')
                    except Exception as e:
                        print(item.qualifiers)
                        raise e


        # Genbank
        outfile = join(snakemake.output.DIR, f"{record.id}.gbk")
        with open(outfile, 'w') as out:
            SeqIO.write(record, out, 'genbank')
    

# Save GC content
df = pd.DataFrame(rows_gc, columns=["NUCCORE_ACC", "NUCCORE_GC", "NUCCORE_Topology"])
df.to_csv(snakemake.output.gc_tab, index=False)

# Save Qualifiers conts
df = pd.DataFrame(rows_qualifiers, 
                  columns=['NUCCORE_ACC']+qualifiers_name+['start', 'end', 'strand'])
df = df.dropna(subset=['protein_id'])
df.to_csv(snakemake.output.proteins_tab, index=False)