#use biopythpom to extract the feature dna sequences foerm genbank
genbank_dir = "/home/davide/Desktop/genomiChro/genbanks_prokka"
nonintergenihce_dir="/home/davide/Desktop/genomiChro/sequenzegeniche"
from Bio import SeqIO
import os

def extract_genic_sequences(genbank_dir, nonintergenihce_dir):
    for file in os.listdir(genbank_dir):
        if file.endswith(".gbk"):
            gen=""
            genoma=SeqIO.parse(os.path.join(genbank_dir, file), 'genbank')
            for seq in genoma:
                for record in seq.features:
                    if record.type in ('CDS','tRNA','rRNA'):
                        gen+=record.extract(seq.seq)
            gen = SeqIO.SeqRecord(gen, id=file[:-4])
            SeqIO.write(gen,os.path.join(nonintergenihce_dir, file[:-4]+"_geniche.fasta"), "fasta")

# Use the function
#extract_genic_sequences("/home/davide/Desktop/genomiChro/genbanks_prokka", "/home/davide/Desktop/genomiChro/sequenzegeniche")
from cercapattern import cercamotivi
motivi = ['.AGT.{5,10}ACT.','TAGT.{7,11}CTA']
motivi_dir="/home/davide/Desktop/genomiChro/sequenzegeniche/motivi"
for file in os.listdir(nonintergenihce_dir):
    if file.endswith(".fasta"):
        out_file=os.path.join(motivi_dir,file[:-6]+".txt")
        cercamotivi(os.path.join(nonintergenihce_dir,file),motivi,out_file)