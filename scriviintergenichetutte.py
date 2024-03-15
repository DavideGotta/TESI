import os
from Bio import SeqIO

from intergeniche import sequenze_intergenichetutte

input_dir= "/home/davide/Desktop/genomiChro/annotati_Refseq"
output_dir = "/home/davide/Desktop/genomiChro/intergeniche_RefSeq"

for file in os.listdir(input_dir):
    input_file = os.path.join(input_dir, file)
    output_file = os.path.join(output_dir, file[:-5]+"_intergen.fasta")
    sequenze_intergenichetutte(input_file, output_file)
    print(f"Found intergenic sequences in {input_file} and saved results to {output_file}")