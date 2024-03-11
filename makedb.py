import subprocess
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
def makeblastdb(input_file, dbtype, output_file):
    """Create a BLAST database from a FASTA file."""
    # Create a BLAST database from a FASTA file
    subprocess.run(['makeblastdb', '-in', input_file, '-dbtype', dbtype, '-out', output_file])
    print(f"Created BLAST database {output_file} from {input_file}")
database_dir= "/home/davide/Downloads/genomiChro/database"
genbank_dir = "/home/davide/Downloads/genomiChro/GenBank"
from estraiproteine import extract_all_proteins_prokka, extract_all_proteins_prokka2

for file in os.listdir(genbank_dir):
    if file.endswith(".gbk"):
        input_file = os.path.join(genbank_dir, file)
        output_file = os.path.join(database_dir, file[:-4] + ".fasta")
        print(output_file)
        extract_all_proteins_prokka(input_file, output_file)
        makeblastdb(output_file, "prot", output_file)