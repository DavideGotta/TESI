from extractproteins import extract_proteins_from_genes
from Bio import SeqIO
import os
gene="lexA"
genbanks_dir = "/home/davide/Desktop/genomiChro/genbanks_prokka"
output_file = "/home/davide/Desktop/genomiChro/lexA_proteins.fa"
for file in os.listdir(genbanks_dir):
    if file.endswith(".gbk"):
        input_file = os.path.join(genbanks_dir, file)
        proteins = extract_proteins_from_genes(input_file, [gene])
        print(f"Extracted {len(proteins)} proteins from {input_file}")
        with open(output_file, "a") as output:
            for protein in proteins:
                SeqIO.write(protein, output, "fasta")