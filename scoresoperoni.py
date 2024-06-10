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
    parser.add_argument('-d', '--intergen_dir', type=str, required=True,
                        help='La directory con le sequenze intergeniche in formato fasta')
    parser.add_argument('-g', '--genome_ref', type=str, required=False,
                        help='Il genoma di interesse su cui fare le analisi, deve essere una sottostringa dei nomi dei '
                             'file nella directory intergeniche(default: CCMEE_29)',
                        default="CCMEE_29")

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
    return sum(1 for i in range(len(s1)) if s1[i] != s2[i])




def refined_score(motivo: Bio.motifs, pid: str, intergenica: str, entropy: list[float], intergeniche_ortologhi:dict, pwm) -> tuple:
    """
    Aggiunge allo score SM qualcosa in base ai motivi che si trovano davanti a geni ortologhi
    :param motivo:  il motivo in formato Bio.motifs
    :param pid:     l'identificatore del gene(Protein ID)
    :param intergenica: la sequenza intergenica upstream a pid nel genoma originale(CCMEE29)
    :param entropy: l'entropia relativa del motivo
    :param intergeniche_ortologhi:  le sequenze intergeniche upstream degli ortologhi del gene
    :param pwm:     la matrice di probabilità delle basi per ogni posizione del motivo
    :return:    lo score della sequenza rispetto al motivo
    """
    info=""
    info_geni=""
    operone= pid.split("|")
    motivo_len = motivo.length
    sm_score = sm(motivo, intergenica, entropy, pwm)
    s = str(sm_score[1])  # sequenza  nell'intergenica di CCMEE29 con maggior score
    max=-100
    for gene in operone:
        if gene not in intergeniche_ortologhi.keys():
            continue
        score_ortologhi = []
        info_geni+=f"Gene {gene}. Score ortologhi= "
        for ceppo,seqs in intergeniche_ortologhi[gene].items():
            info_geni+=f"{ceppo}: {{"
            max_score=-100
            for seq in seqs:
                if len(seq) < motivo_len or "N" in seq:
                    continue
                sm_ortologo = sm(motivo, seq, entropy, pwm)
                score=((motivo_len - Hamming(s, sm_ortologo[1])) / motivo_len )* sm_ortologo[0]
                info_geni+=f"({sm_ortologo[1]},{sm_ortologo[0]:.2f},{score:.2f}), "
                if score > max_score:
                    max_score = score
            info_geni+="}, "
            if max_score==-100:
                max_score=0
            score_ortologhi.append(max_score)
        avg = sum(score_ortologhi) / len(score_ortologhi) if score_ortologhi else 0
        score_ortologhi = [round(x, 2) for x in score_ortologhi]
        info_geni+=f"Max scores degli ortologhi={score_ortologhi}."
        info_geni+=f" Score medio degli ortologhi di {avg:.2f}. "
        if avg>max:
            max=avg
    if max==-100:
        max=0
    return round(sm_score[0],2), round(sm_score[0] + max,2), s, sm_score[2],info_geni


def main():
    args = parse_arguments()
    file_meme_xml = args.meme_xml
    motif_index = args.motif_index
    output_file = args.output_file
    intergen_dir = Path(args.intergen_dir)
    genome_ref = args.genome_ref
    command_line = ET.parse(file_meme_xml).find("model").find("command_line").text
    meme_record = motifs.parse(open(file_meme_xml), 'meme')  # secondo modo per leggere file MEME
    motivo = meme_record[motif_index]
    motivo.pseudocounts = 1
    for intergen in intergen_dir.iterdir():
        if genome_ref in intergen.name:
            intergen = str(intergen)
            break
    else:
        raise FileNotFoundError(f"File {genome_ref} not found in directory {intergen_dir}")
    motivo.background = calculate_background(intergen)
    pwm = motivo.pwm
    data = []
    data.append(command_line)
    entropy = relative_entropy(motivo)
    with open("/home/davide/PycharmProjects/TESI2/dizop.pickle", "rb") as f:
        diz = pickle.load(f)
    start_time = time.time()
    for record in SeqIO.parse(intergen, "fasta"):
        if len(record.seq) < len(motivo) or "N" in record.seq:
            continue
        start = record.description.find("pidCCMEE29=")
        end = record.description.find(",", start)
        pid = record.description[start:end]
        pid=pid.replace("pidCCMEE29=","")

        data.append(
            str(record.id) + "\t" + pid+ "\t" + str(len(record.seq)) + "\t" + "\t".join(map(str, refined_score(motivo, pid, str(record.seq), entropy, diz, pwm))))


    with open(output_file, "w") as f:
        f.writelines("\n".join(data))


if __name__ == "__main__":
    main()
