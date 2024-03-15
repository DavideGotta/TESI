
genomi_dir = "/home/davide/Desktop/genomiChro/annotati_Refseq"
file="Chroococcidiopsis_sp._CCMEE_29_(cyanobacteria)_GCF_023558375.gbff"
from pathlib import Path
import os
file=os.path.join(genomi_dir, file)
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#parse gbff with Biopython and extract fasta with gene in genes
from pprint import pprint
for file in os.listdir(genomi_dir):
    path=os.path.join(genomi_dir, file)
    genoma=SeqIO.parse(path,"genbank")
    i,j=0,0
    go_proc=[]
    for seq in genoma:
        for feature in seq.features:
            if feature.type == "CDS":
                try:
                    #print(feature.qualifiers["gene"])
                    GO=[key for key in feature.qualifiers.keys() if key.startswith("GO")]
                    if "GO_process" in GO:
                        go_proc.append(feature.qualifiers["GO_process"])
                    if not GO:
                       # print("No GO terms found")
                        i+=1
                    else:
                        #print(GO)
                        j+=1
                except KeyError:
                    pass
    print(f"Genes without GO terms: {i}")
    print(f"Genes with GO terms: {j}")
    from collections import Counter
    go_proc = [tuple(lst) for lst in go_proc]
    go_proc = Counter(go_proc)
    pprint(go_proc.most_common(10))
    #calculate p-value of gene set enrichment analysis

