from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
def extract_all_proteins_prokka(input,output):
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
                # Add the protein sequence to the list as a SeqRecord
                proteins.append(
                    SeqRecord(
                        Seq(protein_seq),
                        id=feature.qualifiers["locus_tag"][0],
                        description=feature.qualifiers["product"][0]+" "+feature.qualifiers.get("gene",[' '])[0]+" "+str(feature.location)
                    )
                )

    # Write the protein sequences to a FASTA file
    SeqIO.write(proteins, output, "fasta")
def extract_all_proteins_prokka2(input,output):
    # Create a list to hold the protein sequences
    proteins = []
    genoma = SeqIO.parse(input, "genbank")

    # Loop through the features in the GenBank file
    for feature in genoma.features:
        # Check if the feature is a CDS
        if feature.type == "CDS":
            print(feature)
            # Get the protein sequence
            protein_seq = feature.qualifiers["translation"][0]
            # Add the protein sequence to the list as a SeqRecord
            proteins.append(
                SeqRecord(
                    Seq(protein_seq),
                    id=feature.qualifiers["locus_tag"][0],
                    description=feature.qualifiers["product"][0] + " " + feature.qualifiers.get("gene", [' '])[
                        0] + " " + str(feature.location)
                )
            )

    # Write the protein sequences to a FASTA file
    SeqIO.write(proteins, open(output, "w"), "fasta")
input = "/home/davide/Downloads/genomiChro/GenBank/Chroococcidiopsis_sp._CCMEE_29_(cyanobacteria)_GCF_023558375.1.gbk"
extract_all_proteins_prokka(input,"/home/davide/Downloads/genomiChro/ciao.fasta")
