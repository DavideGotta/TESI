

with open("/home/davide/Downloads/genomiChro/intergeniche_tutte/motivitrovati/Chroococcidiopsis_sp._CCMEE_29_(cyanobacteria)_GCF_023558375.1_intergen_motifs.txt", 'r') as file:
    content = file.read()

geni = []
lines=content.split("\n")
for i, line in enumerate(lines):
    if line.startswith("Risultati per il gene:") and lines[i+2].startswith('Trovato'):
        geni.append(lines[i].split(" ")[4])
print(geni)
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# Open and parse the GenBank file
with open("/home/davide/Downloads/genomiChro/genbanks_prokka/Chroococcidiopsis_sp._CCMEE_29_(cyanobacteria)_GCF_023558375.1.gbk", "r") as file:
    records = SeqIO.parse(file, "genbank")
    proteins = []

    # Iterate over each record in the GenBank file
    for record in records:
        # Iterate over each feature in the record
        for feature in record.features:
            # If the feature is a CDS and its locus_tag is in geni
            if feature.type == "CDS" and feature.qualifiers["locus_tag"][0] in geni:
                # Create a SeqRecord for the protein sequence and add it to the list
                protein_record = SeqRecord(
                    Seq(feature.qualifiers["translation"][0]),
                    id=feature.qualifiers["locus_tag"][0],
                    description=feature.qualifiers["product"][0]+" "+"gene:"+feature.qualifiers.get("gene", ["None"])[0]
                )
                proteins.append(protein_record)

# Write the protein sequences to the output file
SeqIO.write(proteins, "/home/davide/Downloads/genomiChro/genbanks_prokka/CCMEE_proteins_withmotifs.fasta", "fasta")
