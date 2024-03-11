import os
from intergeniche import sequenze_intergenichetutte

input_dir= "/home/davide/Downloads/genomiChro/genbanks_prokka"
output_dir = "/home/davide/Downloads/genomiChro/intergeniche_tutte"

for file in os.listdir(input_dir):
    if not file.endswith(".gbk"):
        continue
    input_file = os.path.join(input_dir, file)
    output_file = os.path.join(output_dir, file[:-4]+"_intergen.fasta")
    sequenze_intergenichetutte(input_file, output_file)
    print(f"Found intergenic sequences in {input_file} and saved results to {output_file}")