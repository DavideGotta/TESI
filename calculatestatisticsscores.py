import argparse
import xml.etree.ElementTree as ET
from collections import Counter
from math import log
from pathlib import Path
import time
import pickle

import Bio
from Bio import SeqIO
from Bio import motifs
import random

def parse_arguments():

    parser = argparse.ArgumentParser(
        description='Dato un motivo crea un file coi motivi trovati nelle intergeniche del genoma di interesse col '
                    'maggiore score SM e quello raffinato rispetto agli ortologhi\nL\'output sarà formato da 5 colonne('
                    'locus_tag,score,score raffinato,motivo trovato, posizione rispetto a TSS) con tante righe quante '
                    'sono le sequenze intergeniche\nEsempio: python3 scores.py -m meme.xml -i 0 -o output.txt -d '
                    'intergeniche -g CCMEE_29',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-m', '--meme_xml', type=str, required=True, help='Il file MEME XML contenente il motivo')
    parser.add_argument('-i', '--motif_index', type=int, required=True, help='L\'indice del motivo nel file MEME XML')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Il file di output dove scrivere i risultati')


    return parser.parse_args()


def calculate_background(file)  -> dict[str, float]:
    """
    Calcola la frequenza delle basi in tutte le sequenze intergeniche  date nel file.fasta
    :param file:    il file fasta con le sequenze intergeniche
    :return:    la frequenza delle basi in tutte le sequenze intergeniche
    """
    sequences = SeqIO.parse(file, "fasta")
    all_sequences = "".join(str(record.seq) for record in sequences)
    nucleotide_counts = Counter(all_sequences)
    total_nucleotides = sum(nucleotide_counts.values())
    nucleotide_frequency = dict()
    for nucleotide, count in nucleotide_counts.items():
        nucleotide_frequency[nucleotide] = count / total_nucleotides
    return nucleotide_frequency


def relative_entropy(motivo: Bio.motifs) -> list[float]:
    """
    Calcola l'entropia relativa per ogni posizione del motivo
    :param motivo:  il motivo in formato Bio.motifs
    :return:    la lista delle entropie relative per ogni posizione del motivo
    """
    q = motivo.background  # frequenze delle basi in tutte le sequenze intergeniche
    pwm = motivo.pwm  # matrice di probabilità delle basi per ogni posizione del motivo
    n = motivo.num_occurrences  # numero di sequenze con cui è stato costruito il motivo
    a = (n + 1) / (n + 4) * log(n + 1) - log(n + 4) - 1 / (n + 4) * sum(log(q[b]) for b in "ACGT") - n / (
            n + 4) * log(min(q.values()))
    entropy = []
    for i in range(motivo.length):
        entropy.append(sum(pwm[b, i] * log(pwm[b, i] / q[b]) for b in "ACGT") / a)
    return entropy


def sm(motivo: Bio.motifs, seq: str, rel_entropy: list, pwm) -> tuple:
    """
    Calcola lo score SM di una sequenza rispetto ad un motivo :param motivo:  il motivo in formato Bio.motifs :param
    seq:     la sequenza su cui calcolare lo score :return:    lo score migliore nella sequenza intergenica rispetto
    al motivo, la sottosequenza a cui corrisponde e la usa poszione rispetto a inizio trascrizione
    """
    q = motivo.background  # frequenze delle basi in tutte le sequenze intergeniche
    max = -float("inf")
    max_i = 0
    for i in range(len(seq) - motivo.length + 1):
        h = seq[i:i + motivo.length]  # sottosequenza di lunghezza del motivo(l-mero della sequenza intergenica)
        score = sum(rel_entropy[i] * log(pwm[h[i], i] / q[h[i]]) for i in range(len(h)))  # score della sottosequenza
        if score > max:
            max = score
            max_i = i
    return max, seq[max_i:max_i + motivo.length], max_i - len(seq)

def Hamming(s1: str, s2: str) -> int:
    """Calcola la distanza di Hamming tra due sequenze"""
    if len(s1) != len(s2):
        raise ValueError("Le sequenze devono avere la stessa lunghezza")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def extract_300_random(lunghezza,cdsCCMEE29,num_seqs=300):
    seqs=[]
    for key in cdsCCMEE29:
        if len(cdsCCMEE29[key])>lunghezza:
            seqs.append((key,cdsCCMEE29[key]))
    random_seqs=[]
    for i in range(num_seqs):
        seq=random.choice(seqs)
        start=random.randint(0,len(seq[1])-lunghezza)
        random_seqs.append((seq[0],seq[1][start:start+lunghezza]))
    return random_seqs

def refined_score_coding(motivo: Bio.motifs, cds: dict, cdsCCMEE29: dict, entropy: list[float], pwm, seq: tuple) -> tuple:
    motivo_len = motivo.length
    score_ortologhi = []
    sm_score = sm(motivo, seq[1], entropy, pwm)
    s = str(sm_score[1])
    if seq[0] in cds:
        for orto in cds[seq[0]]:
            sm_ortologo = sm(motivo, orto, entropy, pwm)
            score=((motivo_len - Hamming(s, sm_ortologo[1])) / motivo_len )* sm_ortologo[0]
            score_ortologhi.append(score)
    avg = sum(score_ortologhi) / len(score_ortologhi) if score_ortologhi else 0

    return sm_score[0], sm_score[0]+avg


def main():
    args = parse_arguments()
    file_meme_xml = args.meme_xml
    motif_index = args.motif_index
    output_file = args.output_file
    meme_record = motifs.parse(open(file_meme_xml), 'meme')  # secondo modo per leggere file MEME
    motivo = meme_record[motif_index]
    motivo.pseudocounts = 1
    file = "/home/davide/Desktop/genomiChro/annotati_Refseq/intergenicheoperoni/CCMEE_29.fasta"
    motivo.background = calculate_background(file)
    pwm = motivo.pwm
    entropy = relative_entropy(motivo)
    with open("/home/davide/PycharmProjects/TESI2/cds.pkl", "rb") as f:
        cds = pickle.load(f)
    with open("/home/davide/PycharmProjects/TESI2/cdsCCMEE29.pkl", "rb") as f:
        cdsCCMEE29 = pickle.load(f)
    seqs = []
    i=0
    # for record in SeqIO.parse(file, "fasta"):
    #     if len(record.seq) < len(motivo) or "N" in record.seq:
    #         continue
    #
    #     i+=1
    #     random_seqs = extract_300_random(len(record.seq), cdsCCMEE29)
    #     seqs.extend(random_seqs)
    scores=[]
    with open(output_file, "w") as f:
        for seq in cdsCCMEE29.items():
            score=refined_score_coding(motivo, cds, cdsCCMEE29, entropy, pwm, seq)
            scores.append(score)
            f.write(f"{score[0]}\t{score[1]}\n")




if __name__ == "__main__":
    main()
