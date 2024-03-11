from extractbesthits import *
import os
dir_genbank = "/home/davide/Downloads/genomiChro/genbanks_prokka"
dir= "/home/davide/Downloads/genomiChro/blast_results_all"
best_dir = "/home/davide/Downloads/genomiChro/blast_results_all/best_hits"
intergeniche_dir = "/home/davide/Downloads/genomiChro/intergeniche_tutte"

for file in os.listdir(dir):
    if file.endswith(".txt"):
        output_file = os.path.join(best_dir, file[:-10] + "_best_hits.txt")
        extract_best_hits(os.path.join(dir, file), output_file)
        print(f"Extracted best hits to {output_file}")

import re
from intergeniche import sequenze_intergeniche
for file in os.listdir(best_dir):

    if file.endswith(".txt"):
        gbk_file = os.path.join(dir_genbank, file[11:-14]+".gbk")
        print(gbk_file)
        with open(os.path.join(best_dir,file), 'r') as f:
            content = f.read()
        query = []
        for i,line in enumerate(content.split("\n")):
            if line.startswith("Query=") and content.split("\n")[i+2].startswith(">"):
                query.append([' '.join(line.split(" ")[1:])+" "+content.split("\n")[i+1]])
        #matches = re.findall(r'>(\S+)', content)
        matches=[line.split(" ")[0][1:] for line in content.split("\n") if line.startswith(">")]

        # Use re.search() to find the expect score and percentage identities
        e_score = re.findall(r"Expect = ([\de-]+)", content)

        perc_id = re.findall(r"Identities = \d+/\d+ \((\d+)%\)", content)
        print(len(matches), len(e_score), len(perc_id), len(query))
        min_len = min(len(matches), len(e_score), len(perc_id), len(query))
        lista = [(matches[i], e_score[i], perc_id[i], query[i]) for i in range(min_len)]
        output_file = os.path.join(intergeniche_dir, file[:-4] + "intergeniche_best_hits.fasta")
        sequenze_intergeniche(gbk_file, lista, output_file)
        print(f"Extracted intergenic sequences to {output_file}")