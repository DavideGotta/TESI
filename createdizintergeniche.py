
from Bio import SeqIO
from pathlib import Path
import pickle
import pandas as pd
# Initialize an empty dictionary
pid_dict = {}
pidsortologhi=pd.read_pickle("ortologhimmseqs.pkl")
#list of the pids fo column CCMEE_29
pids_CCMEE_29=pidsortologhi["CCMEE_29"].values.tolist()
# Directory containing the fasta files
fasta_dir = Path("/home/davide/Desktop/genomiChro/provafgenesboperoni/intergeniche_operoni")

# Iterate over each pid
for pid in pids_CCMEE_29:
    # Iterate over all the fasta files in the directory
    for fasta_file in fasta_dir.iterdir():
        if fasta_file.stem == "CCMEE_29":
            continue
        # Parse the fasta file
        records = SeqIO.parse(str(fasta_file), "fasta")
        # Iterate over each record
        for record in records:
            # Check if the pid is in the description of the record
            if pid in record.description:
                # Add it to the dictionary
                pid_dict.setdefault(pid, {})[fasta_file.stem] = [str(record.seq)]
#save the dictionary to a pickle file
with open("piddict.pkl", "wb") as f:
    pickle.dump(pid_dict, f)