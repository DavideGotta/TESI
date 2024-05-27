import argparse
import xml.etree.ElementTree as ET
from collections import Counter
from math import log
from pathlib import Path
import time

import Bio
from Bio import SeqIO
from Bio import motifs
from Bio.motifs import meme

# Create the parser
parser = argparse.ArgumentParser(
    description='Dato un motivo crea un file coi motivi trovati nelle intergeniche del genoma di interesse col maggiore score SM e quello raffinato rispetto agli ortologhi\nEsempio: python3 scores.py -m meme.xml -i 0 -o output.txt',
    formatter_class=argparse.RawTextHelpFormatter)

# Add the arguments
parser.add_argument('-m', '--meme_xml', type=str, required=True, help='Il file MEME XML contenente il motivo')
parser.add_argument('-i', '--motif_index', type=int, required=True, help='L\'indice del motivo nel file MEME XML')
parser.add_argument('-o', '--output_file', type=str, required=True, help='Il file di output dove scrivere i risultati')
parser.add_argument('-d', '--intergen_dir', type=str, required=True, help='La directory con le sequenze intergeniche in formato fasta')
parser.add_argument('-g', '--genome_ref', type=str, required=False, help='Il genoma di interesse su cui fare le analisi, deve essere una sottostringa dei nomi dei file nella directory intergeniche(default: CCMEE_29)', default="CCMEE_29")

# Parse the arguments
args = parser.parse_args()

# Use the arguments
file_meme_xml = args.meme_xml
motif_index = args.motif_index
output_file = args.output_file
intergen_dir = Path(args.intergen_dir)
genome_ref = args.genome_ref

# file_meme_xml="/home/davide/Desktop/genomiChro/MEME/motivo8_oops/meme.xml"
# meme_record = meme.read(open(file_meme_xml)) # primo modo per leggere file MEM
command_line = ET.parse(file_meme_xml).find("model").find("command_line").text
meme_record = motifs.parse(open(file_meme_xml), 'meme')  # secondo modo per leggere file MEME
motivo = meme_record[motif_index]
motivo.pseudocounts = 1


def calculate_background(file):
    """
    Calcola la frequenza delle basi in tutte le sequenze intergeniche  date nel file.fasta
    :param file:    il file fasta con le sequenze intergeniche
    :return:    la frequenza delle basi in tutte le sequenze intergeniche
    """
    sequences = SeqIO.parse(file, "fasta")
    all_sequences = "".join(str(record.seq) for record in sequences)
    nucleotide_frequency = Counter(all_sequences)
    total_nucleotides = sum(nucleotide_frequency.values())
    for nucleotide, count in nucleotide_frequency.items():
        nucleotide_frequency[nucleotide] = count / total_nucleotides
    return nucleotide_frequency


for intergen in intergen_dir.iterdir():
    if genome_ref in intergen.name:
        intergen = str(intergen)
        break
motivo.background = calculate_background(intergen)


def sm(motivo: Bio.motifs, seq: str):
    """
    Calcola lo score SM di una sequenza rispetto ad un motivo
    :param motivo:  il motivo in formato Bio.motifs
    :param seq:     la sequenza su cui calcolare lo score
    :return:    lo score migliore nella sequenza intergenica rispetto al motivo, la sottosequenza a cui corrisponde e la usa poszione rispetto a inizio trascrizione
    """
    q = motivo.background  # frequenze delle basi in tutte le sequenze intergeniche
    pwm = motivo.pwm  # matrice di probabilità delle basi per ogni posizione del motivo

    def a():
        """
        Calcola il fattore di normalizzazione a
        :return:    il fattore di normalizzazione a
        """
        n = motivo.num_occurrences  # numero di sequenze con cui è stato costruito il motivo
        a = (n + 1) / (n + 4) * log(n + 1) - log(n + 4) - 1 / (n + 4) * sum(log(q[b]) for b in "ACGT") - n / (
                    n + 4) * log(min(q.values()))
        return a

    a = a()  # fattore di normalizzazione a

    def Info(i):
        """ Calcola l'entropia relativa per la posizione i del motivo"""
        somma = sum(pwm[b, i] * log(pwm[b, i] / q[b]) for b in "ACGT")
        return somma / a

    max = -float("inf")
    max_i = 0
    for i in range(len(seq) - motivo.length + 1):
        h = seq[i:i + motivo.length]  # sottosequenza di lunghezza del motivo(l-mero della sequenza intergenica)
        score = sum(Info(i) * log(pwm[h[i], i] / q[h[i]]) for i in range(len(h)))  # score della sottosequenza
        if score > max:
            max = score
            max_i = i
    return max, seq[max_i:max_i + motivo.length], max_i - len(seq)


# @title Funzione che calcola lo score raffinato in base agli ortologhi


def refined_score(motivo: Bio.motifs, pid: str, intergenica):
    """
    Aggiunge allo score SM qualcosa in base ai motivi che si trovano davanti a geni ortologhi
    :param motivo:  il motivo in formato Bio.motifs
    :param pid:     l'identificatore del gene(Protein ID)
    :param intergenica: la sequenza intergenica upstream a pid nel genoma originale(CCMEE29)
    :return:    lo score della sequenza rispetto al motivo
    """

    def Hamming(s1: str, s2: str):
        """Calcola la distanza di Hamming tra due sequenze"""

        return sum(1 for i in range(len(s1)) if s1[i] != s2[i])

    def ortologhi(pid: str, intergen_dir: Path):
        """
        Restituisce la lista delle sequenze intergeniche upstream degli ortologhi del gene
        :param pid: l'identificatore del gene(Protein ID)
        :return:    la lista delle sequenze intergeniche upstream degli ortologhi del gene
        """
        seqs = []
        for file in intergen_dir.iterdir():
            if genome_ref not in file.name:
                for record in SeqIO.parse(file, "fasta"):
                    if pid in record.description:
                        seqs.append(record.seq)
        return seqs

    intergeniche_ortologhi = ortologhi(pid, intergen_dir)
    motivo_len = motivo.length
    sm_score = sm(motivo, intergenica)
    s = str(sm_score[1])  # sequenza  nell'intergenica di CCMEE29 con maggior score

    score_ortologhi = [(motivo_len - Hamming(s, sm(motivo, seq)[1])) / motivo_len * sm(motivo, seq)[0] for seq in
                       intergeniche_ortologhi if len(seq) >= motivo_len and "N" not in seq]
    avg = sum(score_ortologhi) / len(score_ortologhi) if len(score_ortologhi) > 0 else 0
    return sm_score[0], sm_score[0] + avg, s, sm_score[2]

def extract_coding_seq(genoma: str):
    """
    Estrae tutte le sequenze codificanti dal genoma genbank
    """
    for record in SeqIO.parse(genoma, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                yield (feature.location.extract(record.seq), feature.qualifiers["locus_tag"][0])


def countseqs(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
def LOR(score):

    num = sum(1 for record in SeqIO.parse(intergen, "fasta") if sm(motivo, record.seq)[0] > score) / countseqs(intergen)
    coding = extract_coding_seq(
        "/home/davide/Desktop/genomiChro/annotati_Refseq/Chroococcidiopsis_sp._CCMEE_29_(cyanobacteria)_GCF_023558375.gbff")
    den = sum(1 for seq in coding if sm(motivo, seq[0]) > score)
    if den == 0:
        return 0
    else:
        return log(num / den)

data= []

data.append(command_line)
for record in SeqIO.parse(intergen, "fasta"):
    if len(record.seq) < len(motivo) or "N" in record.seq:
        continue
    start = record.description.find("WP_")
    end = record.description.find("'", start)
    pid = record.description[start:end]
    data.append(str(record.id) + "\t" + "\t".join(map(str, refined_score(motivo, pid, record.seq))))
start_time = time.time()
with open(output_file, "w") as f:
    f.writelines("\n".join(data))
end_time = time.time()
print(f"Execution time: {end_time - start_time} seconds")
