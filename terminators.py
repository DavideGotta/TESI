from Bio import SeqIO

gff_file = "/home/davide/Desktop/genomiChro/gffs/refseq/genomic.gff"
import pandas as pd
df=pd.read_csv(gff_file, sep="\t", comment="#", header=None)
df.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
print(df.head())