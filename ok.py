import os
import sys
from Bio import SeqIO
geni=['lexA','recA','uvrB','ssb']
dir="/home/davide/Downloads/genomiChro/intergeniche_best_hits"
#open all fasta files in the directory
for file in os.listdir(dir):
    with open("prova.fasta", "a") as f:
        if file.endswith(".fasta"):
            #use biopython to read the fasta file
            for record in SeqIO.parse(os.path.join(dir, file), "fasta"):
                #if the record description has one of the genes of interest, write the record to a new file fasta
                for gene in geni:
                    if gene in record.description:
                        SeqIO.write(record, f, "fasta")