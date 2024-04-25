
import subprocess
import pathlib
from pathlib import Path
dir=Path("/home/davide/Desktop/genomiChro/proteine_RefSeq")
db_name=Path("/home/davide/Desktop/genomiChro/databases/NostocPCC7120")
output_dir=Path("/home/davide/Desktop/genomiChro/blastvsPCC7120")
#iterate dir
for fasta in dir.glob("*.faa"):

    #how to access name of file without all the path

    output_file=output_dir/f"{fasta.stem}Vs{db_name.stem}.txt"

    subprocess.run(["blastp", "-query", fasta, "-db", db_name, "-out", output_file, "-outfmt", "6"  ])
    print(f"BLAST vs {db_name} completed for {fasta.stem}")