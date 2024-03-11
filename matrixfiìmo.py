"""
CREAZUIONE DI UNA MATRICE IN CUI LE RIGHE SONO I GENI E LE COLONNE SONO I GENOMI E I VALORI SONO I MOTIVI TROVATI(E LA POSIZIONE) CON FIMO DEL MOTIVO DI 16BP TROVATO CON MEME

CREAZIONE DI GRAFICO A TORTA CHE MOSTRA LA FREQUENZA DEI CODONI DI START(IN TEORIA AVEVO LETTO CHE PER GENI DI RIPARO SPESSO SI HA UN CODONE DI INIZIO ALTERNATIVO)
"""

# Import necessary libraries
import os
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
codoni=dict()
# Function to get the length of the intergenic sequence
def get_intergenic_length(fimo_file, locus_tag):
    """
    This function returns the length of the intergenic sequence for a given locus tag.

    Parameters:
    fimo_file (str): The name of the FIMO file.
    locus_tag (str): The locus tag of the gene.

    Returns:
    int: The length of the intergenic sequence.
    """
    intergen_dir = "/home/davide/Desktop/genomiChro/intergeniche_tutte"
    for fasta in os.listdir(intergen_dir):
        if fimo_file[4:-10] in fasta:
            for rec in SeqIO.parse(os.path.join(intergen_dir,fasta), 'fasta'):
                if locus_tag in rec.description:
                    return len(rec.seq)
dir="/home/davide/Desktop/genomiChro/fimo"
genbank_dir = "/home/davide/Desktop/genomiChro/genbanks_prokka"
# Initialize DataFrame
dataframe= pd.DataFrame(columns=['gene'])

# Loop through all files in the directory
for file in os.listdir(dir):
    if file.startswith("."):  # Skip hidden files
        continue
    # Read data from file
    data = pd.read_csv(os.path.join(dir, file), delimiter='\t', index_col=None)
    data = data[:-4]  # Delete last 4 rows that contain FIMO algorithm specs
    # Reorder columns with sequence_name first
    data = data[['sequence_name', 'start','matched_sequence']]
    data['start'] = data['start'].astype(int)
    # Apply function to sequence_name column
    lunghezze_intergeniche = data['sequence_name'].apply(lambda x: get_intergenic_length(file,x))
    data['start'] =  data['start'] - lunghezze_intergeniche

    # Loop through all files in the genbank directory
    for file2 in os.listdir(genbank_dir):
        if file[4:-10] in file2:
            genoma=SeqIO.parse(os.path.join(genbank_dir, file2), 'genbank')
            for seq in genoma:
                for record in seq.features:
                    if record.type=='CDS' and record.qualifiers['locus_tag'][0] in data['sequence_name'].values:
                        if record.location.strand == 1:
                            codone = seq.seq[record.location.start:record.location.start+3]
                            codoni[codone] = codoni.get(codone, 0) + 1
                        elif record.location.strand == -1 :
                            codone = seq.seq[record.location.end-3:record.location.end].reverse_complement()
                            codoni[codone] = codoni.get(codone, 0) + 1
                        if record.qualifiers.get('gene',0):
                            gene=record.qualifiers['gene'][0]
                            locus_tag=record.qualifiers['locus_tag'][0]
                            data['sequence_name'] = data['sequence_name'].replace(record.qualifiers['locus_tag'][0], gene)
                        elif record.qualifiers['product']!=['hypothetical protein']:
                            gene=record.qualifiers['product'][0]
                            locus_tag=record.qualifiers['locus_tag'][0]
                            data['sequence_name'] = data['sequence_name'].replace(record.qualifiers['locus_tag'][0], gene)
                        else:
                            # Delete the row that contains the sequence name
                            data = data[data.sequence_name != record.qualifiers['locus_tag'][0]]
    data['matched_sequence'] = data['matched_sequence'] + "(" + data['start'].astype(str) + ")"
    data = data.drop(columns=['start'])
    data.rename(columns={'sequence_name':'gene'}, inplace=True)
    for file2 in os.listdir(genbank_dir):
        if file[4:-10] in file2:
            name = file2[file2.find("_")+1:file2.find("(")-1]
    data.rename(columns={'matched_sequence': name}, inplace=True)
    dataframe = pd.merge(dataframe, data, on='gene', how='outer')
    df = pd.DataFrame(list(codoni.items()), columns=['Start Codon', 'Frequency'])
    df = df.sort_values(by='Frequency', ascending=False)
    plt.pie(df['Frequency'], labels=df['Start Codon'], autopct='%1.1f%%')
    plt.title('Start Codon Frequency')
    plt.show()

# Reorder rows sort for number of Nan values
dataframe = dataframe.reindex(dataframe.isnull().sum(axis=1).sort_values().index)

# Save DataFrame to excel
dataframe.to_excel('/home/davide/Downloads/geniconMotivoMEMEinGenomiChro.xlsx', index=False)

# Save DataFrame to csv
dataframe.to_csv('/home/davide/Downloads/geniconMotivoMEMEinGenomiChro.csv', index=False)