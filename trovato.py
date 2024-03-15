motivi_dir=("/home/davide/Desktop/genomiChro/intergeniche_RefSeq/motivo")
genbank_dir = "/home/davide/Desktop/genomiChro/annotati_Refseq"
import re
import os
for file in os.listdir(genbank_dir):
    if file.endswith(".gbff"):
        #search locus_tag="esempio" and extract the gene name
        with open(os.path.join(genbank_dir, file), 'r') as f:
            x = re.search(r'locus_tag="(.+?)"', f.read()).group(1)
            print(f"'{x.split('_')[0]}' : '{file}' ,")