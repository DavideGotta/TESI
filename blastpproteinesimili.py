import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# Specify the gene list file
#fasta = "/home/davide/Downloads/genomiChro/genbanks_prokka/CCMEE_proteins_withmotifs.fasta"
fasta="/home/davide/Desktop/genomiChro/proteineriparo.fasta"
# Specify the directory containing the protein databases
db_dir = "/home/davide/Desktop/genomiChro/databases"

ris_dir = "/home/davide/Desktop/genomiChro/riparo"
# Get a list of all protein databases in the directory
db_files = os.listdir(db_dir)

# Iterate over each protein database in the directory
for db_file in db_files:
    # If the file is a protein database
    if db_file.endswith(".fasta"):
        print(db_file)
        # Construct the path to the protein database
        db_path = os.path.join(db_dir, db_file)
        #remove extension from db_file
        db_file=db_file[:-6]
        # Construct the output file path
        output_file = os.path.join(ris_dir, f"riparo_{db_file}.txt")

        # Execute the blastp command
        subprocess.run(["blastp", "-query", fasta, "-db", db_path, "-out", output_file])