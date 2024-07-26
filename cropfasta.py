from Bio import SeqIO
file="/home/davide/Desktop/genomiChro/fasta/Chroococcidiopsis_sp_CCMEEE_29_chromosome.fa"
#crop fasta to stay under 5 Mb of data
for record in SeqIO.parse(file, "fasta"):
    if len(record.seq) > 4800000:
        record.seq = record.seq[:4800000]
    print(f"crop to 5 Mb: {record.id}")
    with open(f"/home/davide/Desktop/genomiChro/fasta/cropped/{record.id}.fa", "w") as f:
        SeqIO.write(record, f, "fasta")