from cercapattern import cercamotivi
import os

input_dir= "/home/davide/Downloads/genomiChro/intergeniche_tutte/intergeniche_best"
output_dir = "/home/davide/Downloads/genomiChro/intergeniche_tutte/motivibest"
patterns = ['.AGT.{3,11}ACT.', '.TAGT.{3,11}ATC.', '.AG.{3,11}ACT.', '.AGT.{3,11}CT.' ]

for file in os.listdir(input_dir):
    if not file.endswith(".fasta"):
        continue
    input_file = os.path.join(input_dir, file)
    output_file = os.path.join(output_dir, file[:-6]+"_motifs.txt")
    cercapatterns(input_file, patterns, output_file)
    print(f"Found patterns in {input_file} and saved results to {output_file}")



# Specify the directory containing the text files
input_dir = "/home/davide/Downloads/genomiChro/intergeniche_best_hits/motivitrovati"

# Specify the output file
output_file = "/home/davide/Downloads/genomiChro/intergeniche_best_hits/motivitrovati/joined.txt"

# Use the cat command to join all text files into a single file
os.system(f"cat {input_dir}/*.txt > {output_file}")
import re
import collections

# Specify the output file
output_file = "/home/davide/Downloads/genomiChro/intergeniche_best_hits/motivitrovati/joined.txt"

# Open the output file and read its contents
with open(output_file, 'r') as file:
    lines = file.readlines()

# Initialize a list to store the strings
strings = []

# Iterate over each line in the file
for i in range(len(lines)):
    # If the line contains a string inside [] square brackets
    match = re.search(r'\'(.*?)\'', lines[i])
    if match and i + 2 < len(lines) and lines[i + 2].startswith('Trovato'):
        # Add the string to the list
        strings.append(match.group(1))

# Count how many times each unique string occurs
counts = collections.Counter(strings)

#sort for biggest number of occurences
sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
# Print the counts

locus_tags = [string.split(" ")[0] for string, count in sorted_counts]
new_counts = dict(zip(locus_tags,  [count for string, count in sorted_counts]))

from Bio import SeqIO
# List of locus_tags you are interested in
#infile="Chroococcidiopsis_sp._CCMEE_29_(cyanobacteria)_GCF_023558375.1.gbk"
infile="prokkaCCMEE29.gbk"
indir="/home/davide/Downloads/PROKKA_02102024"
#indir="/home/davide/Downloads/genomiChro/genbanks_prokka"

# Open and parse the GenBank file
annotati={}
with open(os.path.join(indir, infile), "r") as file:
    records = SeqIO.parse(file, "genbank")
    proteins = []

    # Iterate over each record in the GenBank file
    for record in records:
        # Iterate over each feature in the record
        for feature in record.features:
            # If the feature is a CDS and its locus_tag is in locus_tags
            if feature.type == "CDS" and feature.qualifiers["locus_tag"][0] in locus_tags:
                annotati[feature.qualifiers["locus_tag"][0]] ="product = "+feature.qualifiers["product"][0]+" "+"gene = "+feature.qualifiers.get("gene", ["None"])[0]+" "+str(feature.location)

# Define the output file path
output_file_path = '/home/davide/Downloads/genomiChro/intergeniche_best_hits/motivitrovati/geniperoccorrenze.txt'
with open(output_file_path, 'w') as outfile:
    for key, value in new_counts.items():
        try:
            outfile.write(f"{key} {annotati[key]} Occorrenze = {value}\n")
        except KeyError:
            print(f"Key {key} not found in the GenBank file")

