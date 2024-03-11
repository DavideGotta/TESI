#statistica dei codoni di inizio
genbank_dir = "/home/davide/Desktop/genomiChro/genbanks_prokka"
from Bio import SeqIO
#count how many times each start codon is found in the genome
from collections import Counter
import os
codoni = dict()
for file in os.listdir(genbank_dir):

    if file.endswith(".gbk"):
        print(file)
        genoma=SeqIO.parse(os.path.join(genbank_dir, file), 'genbank')
        for seq in genoma:
            for record in seq.features:
                if record.type=='CDS':
                    print(record.qualifiers['codon_start']) if record.qualifiers['codon_start']!=['1'] else None
                    if record.location.strand == 1:
                        codon = str(seq.seq[record.location.start:record.location.start+3])
                        codoni[codon] = codoni.get(codon, 0) + 1
                        #print(seq.seq[record.location.start:record.location.start+3])
                    else:
                        codon = str(seq.seq[record.location.end-3:record.location.end].reverse_complement())
                        codoni[codon] = codoni.get(codon, 0) + 1
                        #print(seq.seq[record.location.end-3:record.location.end].reverse_complement())
        import matplotlib.pyplot as plt
        import pandas as pd

        df = pd.DataFrame(list(codoni.items()), columns=['Start Codon', 'Frequency'])
        df = df.sort_values(by='Frequency', ascending=False)
        plt.pie(df['Frequency'], labels=df['Start Codon'], autopct='%1.1f%%')
        plt.title('Start Codon Frequency')
        plt.show()

