
from Bio import SeqIO
def clean_fasta(input,output):
    records=[]
    for record in SeqIO.parse(input,"fasta"):
        print(record.description)
        print(record.id)
        record.id=record.id.split("|")[0]
        record.description=record.id
        records.append(record)
    SeqIO.write(records,output,"fasta")

input="es.fasta"
output="clean.fasta"
clean_fasta(input,output)