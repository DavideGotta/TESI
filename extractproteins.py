#in realtà non serve a niente, basta usare il file .faa che prokka genera
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_all_proteins_prokka(input, output):
    # Create a list to hold the protein sequences
    proteins = []
    genoma = SeqIO.parse(input, "genbank")
    # Loop through the features in the GenBank file
    for seq in genoma:
        for feature in seq.features:
            # Check if the feature is a CDS
            if feature.type == "CDS":
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
    SeqIO.write(proteins, output, "fasta")

def extract_proteins_from_genes(input, genenames):
    # Create a list to hold the protein sequences
    proteins = []
    genoma = SeqIO.parse(input, "genbank")
    # Loop through the features in the GenBank file
    for seq in genoma:
        for feature in seq.features:
            # Check if the feature is a CDS
            if feature.type == "CDS":
                # Get the protein sequence
                protein_seq = feature.qualifiers["translation"][0]
                # Add the protein sequence to the list as a SeqRecord
                if any(gene in feature.qualifiers.get('gene',[''])[0] for gene in genenames):
                    proteins.append(
                        SeqRecord(
                            Seq(protein_seq),
                            id=feature.qualifiers["locus_tag"][0]+" "+input.split("/")[-1][:-4],
                            description=feature.qualifiers["product"][0] + " " + feature.qualifiers.get("gene", [' '])[
                                0] + " " + str(feature.location)
                        )
                    )

    # Write the protein sequences to a FASTA file
    return proteins
