file="/home/davide/Downloads/PCC6803GT.gbk"
from Bio import SeqIO
from Bio.Seq import Seq
#import seqrecord
from Bio.SeqRecord import SeqRecord
genoma = SeqIO.parse(file, "genbank")
seqs = []
for seq in genoma:
    for record in seq.features:
        if record.type == "CDS":
            if record.qualifiers.get("protein_id", [""])[0] != "":
                seqrecord = SeqRecord(
                    Seq(record.qualifiers["translation"][0]),
                    id=record.qualifiers["protein_id"][0],
                    description=record.qualifiers.get("product",[""])[0] + " " + record.qualifiers.get("gene", [''])[0] +" "+ record.qualifiers["locus_tag"][0]+"(" + record.qualifiers.get("note",[""])[0]+") "+ str(record.location)
                )
                seqs.append(seqrecord)
from Bio import SeqIO
SeqIO.write(seqs, "/home/davide/Desktop/genomiChro/proteinePCC6803/PCC6803.fasta", "fasta")