import os
import subprocess

# Specify the gene list file
fasta = "/home/davide/Desktop/genomiChro/genbanks_prokka/CCMEE_proteins_withmotifs.fasta"
#fasta="/home/davide/Desktop/genomiChro/proteineriparo.fasta"
# Specify the directory containing the protein databases
db_dir = "/home/davide/Desktop/genomiChro/databases
# Specify the directory to save the BLAST results
ris_dir = "/home/davide/Desktop/genomiChro/blast_results"

# Iterate over each protein database in the directory
for db_file in os.listdir(db_dir):
    if db_file.endswith(".fasta"):
        db_path = os.path.join(db_dir, db_file)
        db_file=db_file[:-6]
        output_file = os.path.join(ris_dir, f"riparo_{db_file}.txt")
        subprocess.run(["blastp", "-query", fasta, "-db", db_path, "-out", output_file])