import Bio.SeqIO
import pyrodigal

records = [record for record in Bio.SeqIO.parse(snakemake.input.fasta, "fasta")]

orf_finder = pyrodigal.GeneFinder(meta=True)
for i, pred in enumerate(orf_finder.find_genes(bytes(records[0].seq))):
    print(f">{records[0].id}_{i+1}")
    print(pred.translate())