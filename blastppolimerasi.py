import os
import subprocess

# Specify the gene list file
fasta = "proteineAnaCd.fasta"
# Specify the directory to save the BLAST results
ris_dir = "/home/davide/Desktop/genomiChro/"
db_dir = "/home/davide/Desktop/genomiChro/database_RefSeq"
# Iterate over each protein database in the directory
import os
os.chdir(db_dir)
for db_file in os.listdir(db_dir):
    if db_file.endswith(".pin"):
        db_name=db_file[:-4]
        output_file = os.path.join(ris_dir, f"{db_file}.txt")
        subprocess.run(["blastp", "-query", fasta, "-db", db_name, "-out", output_file])