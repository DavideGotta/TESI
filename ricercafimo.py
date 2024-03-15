intergen_dir = "/home/davide/Desktop/genomiChro/intergeniche_RefSeq"
import os
import subprocess
for fasta in os.listdir(intergen_dir):
    input_file = os.path.join(intergen_dir, fasta)
    output_file = os.path.join(intergen_dir, fasta[:-6] + ".blastdb")
    subprocess.run(['makeblastdb', '-in', input_file, '-dbtype', 'nucl', '-out', output_file])
    print(f"Created BLAST database {output_file} from {input_file}")