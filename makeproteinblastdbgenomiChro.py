import os
import subprocess


def makeblastdb(input_file, dbtype, output_file):
    """Create a BLAST database from a FASTA file."""
    subprocess.run(['makeblastdb', '-in', input_file, '-dbtype', dbtype, '-out', output_file])
    print(f"Created BLAST database {output_file} from {input_file}")


database_dir = "/home/davide/Desktop/genomiChro/database_RefSeq"
genbank_dir = "/home/davide/Desktop/genomiChro/proteine_RefSeq"
from extractproteins import extract_all_proteins_prokka
import os
os.chdir(database_dir)
for file in os.listdir(genbank_dir):
    if file.endswith(".faa"):
        input_file = os.path.join(genbank_dir, file)
        output_file = file[:-4]
        makeblastdb(input_file, "prot", output_file)
