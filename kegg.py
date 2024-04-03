from Bio import KEGG
from Bio.KEGG import REST
from Bio import SeqIO
#import SeqRecord
from Bio.SeqRecord import SeqRecord


lista="""alr1742 all4337 all0005 alr0880 alr0803 alr1827 alr1890 all0865  all1864 alr0535 alr0882 alr2425 alr1208 all0579 all7597 all0329 all0258 alr0806 all4391 alr0782 all4688 alr5302  all4894 alr0052 alr2205 alr3829 all5062 all3909 alr4853  all3797 all2110 alr0528 all4050 all0136 all4377 all7618 all7619 alr2493 alr7622 alr4267 all0121"""
proteine=["ana:"+x for x in lista.split()]
seqs=[]
for proteina in proteine:
    #write the protein sequence to a file
    print(f"Downloading {proteina}")
    try:
        proteina=REST.kegg_get(proteina, "aaseq")
        seq = SeqIO.read(proteina, "fasta")
        seqs.append(seq)
    except Exception as e:
        print(e)
        continue
    #read like seqrecord

SeqIO.write(seqs, "proteineAna.fasta", "fasta")
