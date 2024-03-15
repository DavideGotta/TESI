from Bio import SeqIO
import re


def cercamotivi(input_filename, patterns, output_filename):
    results = {}
    diz = {pattern: 0 for pattern in patterns}
    len_tot=0
    with open(input_filename, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            len_tot+=len(record.seq)
            trovati = set()
            results[record.description] = []
            #diz that count how many each pattern is found in all the genome
            for pattern in patterns:
                matches = [(m.start(), m.group()) for m in re.finditer(pattern, str(record.seq))]
                copia = matches.copy()
                diz[pattern] += len(matches)
                for match in matches:
                    if match not in trovati:
                        trovati.add(match)
                    else:
                        copia.remove(match)
                if copia:
                    results[record.description].append((len(record.seq), pattern, copia))

    with open(output_filename, 'w') as outfile:
        for gene, values in results.items():
            outfile.write(f"Risultati per il gene: {gene}\n"+'\n')
            for length, pattern, matches in values:
                outfile.write(f"Trovato/i {len(matches)} occorrenze del motivo {pattern}\n")
                for position, match_str in matches:
                    match_str = match_str[:1].lower() + match_str[1:-1] + match_str[-1:].lower()
                    outfile.write(f"\t-il motivo {match_str} alla posizione {position-length}\n")
            outfile.write('\n'+'*'*100+'\n')
        outfile.write(f"Numero di occorrenze per ogni motivo\n")
        for pattern, count in diz.items():
            outfile.write(f"{pattern} = {count}, normalizzato = {count/len_tot}\n")
        outfile.write(f"Lunghezza totale:{len_tot}")


motivi = ['.AGT.{8}ACT.','TAGT.{9}CTA']
#, '.TAGT.{3,11}ATC.', '.AG.{3,11}ACT.', '.AGT.{3,11}CT.'

import os
# Specify the input directory
dir="/home/davide/Desktop/genomiChro/intergeniche_RefSeq"
# Specify the output directory
output_dir = "/home/davide/Desktop/genomiChro/intergeniche_RefSeq/motivo"
for file in os.listdir(dir):
    #if file is directory or hidden file skip
    if not file.endswith(".fasta"):
        continue
    input_file = os.path.join(dir, file)
    output_file = os.path.join(output_dir, file[:-3]+"_motifs.txt")
    cercamotivi(input_file, motivi, output_file)
    print(f"Found patterns in {input_file} and saved results to {output_file}")
