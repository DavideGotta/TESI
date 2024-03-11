file="/home/davide/geniriparo.txt"
geni=open(file,"r").read().split("\n")
geni=[gene.strip() for gene in geni if gene!=""]
genbank="/home/davide/Desktop/genomiChro/genbanks_prokka/Chroococcidiopsis_sp._CCMEE_29_cyanobacteria_GCF_023558375.1.gbk"
#parse genbank with Biopython and extract fasta with gene in genes
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
def extract_proteins_prokka(input,output,genes):
    # Create a list to hold the protein sequences
    proteins = []
    genoma = SeqIO.parse(input, "genbank")
    # Loop through the features in the GenBank file
    for x in genoma:
        for feature in x.features:
            # Check if the feature is a CDS
            if feature.type == "CDS":
                # Get the protein sequence
                protein_seq = feature.qualifiers["translation"][0]
                gene = feature.qualifiers.get("gene",[''])[0]
                # Add the protein sequence to the list as a SeqRecord
                if gene in genes:
                    proteins.append(
                        SeqRecord(
                            Seq(protein_seq),
                            id=feature.qualifiers["locus_tag"][0],
                            description=feature.qualifiers["product"][0]+" "+gene+" "+str(feature.location)
                        )
                    )
    # Write the protein sequences to a FASTA file
    SeqIO.write(proteins, output, "fasta")
extract_proteins_prokka(genbank,"/home/davide/Desktop/genomiChro/proteineriparo.fasta",geni)