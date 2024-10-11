import argparse
import xml.etree.ElementTree as ET
from collections import Counter
from math import log
import time
import pickle
import random
from multiprocessing import Pool, cpu_count

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
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Il file di output dove scrivere i risultati')
    parser.add_argument('-d', '--intergen_dir', type=str, required=False,
                        help='La directory con le sequenze intergeniche in formato fasta')
    parser.add_argument('-g', '--genome_ref', type=str, required=False,
                        help='Il genoma di interesse su cui fare le analisi, deve essere una sottostringa dei nomi dei '
                             'file nella directory intergeniche(default: CCMEE_29)',
                        default="CCMEE_29")

    return parser.parse_args()


def calculate_background(file) -> dict[str, float]:
    sequences = SeqIO.parse(file, "fasta")
    all_sequences = "".join(str(record.seq) for record in sequences)
    nucleotide_counts = Counter(all_sequences)
    total_nucleotides = sum(nucleotide_counts.values())
    nucleotide_frequency = dict()
    for nucleotide, count in nucleotide_counts.items():
        nucleotide_frequency[nucleotide] = count / total_nucleotides
    return nucleotide_frequency


def relative_entropy(motivo: Bio.motifs) -> list[float]:
    q = motivo.background
    pwm = motivo.pwm
    n = motivo.num_occurrences
    a = (n + 1) / (n + 4) * log(n + 1) - log(n + 4) - 1 / (n + 4) * sum(log(q[b]) for b in "ACGT") - n / (
            n + 4) * log(min(q.values()))
    entropy = []
    for i in range(motivo.length):
        entropy.append(sum(pwm[b, i] * log(pwm[b, i] / q[b]) for b in "ACGT") / a)
    return entropy


def sm(motivo: Bio.motifs, seq: str, rel_entropy: list, pwm) -> tuple:
    q = motivo.background
    max = -float("inf")
    max_i = 0
    for i in range(len(seq) - motivo.length + 1):
        h = seq[i:i + motivo.length]
        score = sum(rel_entropy[i] * log(pwm[h[i], i] / q[h[i]]) for i in range(len(h)))
        if score > max:
            max = score
            max_i = i
    return max, seq[max_i:max_i + motivo.length], max_i - len(seq)


def Hamming(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        raise ValueError("Le sequenze devono avere la stessa lunghezza")
    return sum(1 for i in range(len(s1)) if s1[i] != s2[i])


def refined_score(motivo: Bio.motifs, pid: str, cds: str, entropy: list[float], scores_ortologhi: dict, pwm) -> tuple:
    operone = pid.split("|")
    motivo_len = motivo.length
    sm_score = sm(motivo, cds, entropy, pwm)
    s = str(sm_score[1])
    max = -100
    for gene in operone:
        if gene not in scores_ortologhi.keys():
            continue
        score_ortologhi = []
        for sm_ortologo in scores_ortologhi[gene]:
            score = ((motivo_len - Hamming(s, sm_ortologo[1])) / motivo_len) * sm_ortologo[0]
            score_ortologhi.append(score)
        avg = sum(score_ortologhi) / len(score_ortologhi) if score_ortologhi else 0
        if avg > max:
            max = avg
    if max == -100:
        max = 0
    return round(sm_score[0], 2), round(sm_score[0] + max, 2), s, sm_score[2]


def extract_300_random(lunghezza, cdsCCMEE29, num_seqs=300):
    seqs = []
    for key in cdsCCMEE29:
        if len(cdsCCMEE29[key]) > lunghezza:
            seqs.append((key, cdsCCMEE29[key]))
    random_seqs = []
    for i in range(num_seqs):
        seq = random.choice(seqs)
        start = random.randint(0, len(seq[1]) - lunghezza)
        random_seqs.append(seq[1][start:start + lunghezza])
    return random_seqs


def process_intergenic(info, motivo, entropy, scores, pwm, cdsCCMEE29):
    lunghezza, descrizione = info
    random_seqs = extract_300_random(lunghezza, cdsCCMEE29)
    start = descrizione.find("pids=")
    end = descrizione.find(",", start)
    pids = descrizione[start:end]
    pids = pids.replace("pids=", "")
    results = []
    for seq in random_seqs:
        results.append(pids + "\t" + "\t".join(map(str, refined_score(motivo, pids, seq, entropy, scores, pwm))))
    return results


def main():
    args = parse_arguments()
    file_meme_xml = args.meme_xml
    motif_index = args.motif_index
    output_file = args.output_file
    command_line = ET.parse(file_meme_xml).find("model").find("command_line").text
    meme_record = motifs.parse(open(file_meme_xml), 'meme')
    motivo = meme_record[motif_index]
    motivo.pseudocounts = 1
    intergen = "/home/davide/Desktop/genomiChro/provafgenesboperoni/intergeniche_operoni/CCMEE_29.fasta"
    motivo.background = calculate_background(intergen)
    pwm = motivo.pwm
    data = []
    data.append(command_line)
    with open("scorescoding.pkl", "rb") as f:
        scores = pickle.load(f)
    with open("/home/davide/PycharmProjects/TESI2/cdsCCMEE29.pkl", "rb") as f:
        cdsCCMEE29 = pickle.load(f)
    entropy = relative_entropy(motivo)
    info_intergeniche = [(len(record.seq), record.description) for record in SeqIO.parse(intergen, "fasta")]
    start_time = time.time()
    with Pool(cpu_count()) as pool:
        results = pool.starmap(process_intergenic, [(info, motivo, entropy, scores, pwm, cdsCCMEE29) for info in info_intergeniche])
    end_time = time.time()
    print(f"Elapsed time: {end_time - start_time:.2f} seconds")
    for result in results:
        data.extend(result)

    with open("/home/davide/Documents/datastats.pkl", "wb") as f:
        pickle.dump(data, f)
    with open(output_file, "w") as f:
        f.writelines("\n".join(data))


if __name__ == "__main__":
    main()