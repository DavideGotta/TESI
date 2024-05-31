import argparse
import xml.etree.ElementTree as ET
from collections import Counter
from math import log
from pathlib import Path
import time
import numpy as np

import Bio
from Bio import SeqIO
from Bio import motifs
from Bio.motifs import meme


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Given a motif, create a file with motifs found in the intergenic regions of the genome of interest with the highest SM score and the refined score against orthologs.\n"
            "The output will consist of 5 columns (locus_tag, score, refined score, found motif, position relative to TSS) with as many rows as there are intergenic sequences.\n"
            "Example: python3 scores.py -m meme.xml -i 0 -o output.txt -d intergenic -g CCMEE_29"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument('-m', '--meme_xml', type=str, required=True, help='The MEME XML file containing the motif')
    parser.add_argument('-i', '--motif_index', type=int, required=True, help='The motif index in the MEME XML file')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='The output file to write the results to')
    parser.add_argument('-d', '--intergen_dir', type=str, required=True,
                        help='The directory with intergenic sequences in FASTA format')
    parser.add_argument('-g', '--genome_ref', type=str, required=False,
                        help='The genome of interest to analyze (default: CCMEE_29)', default="CCMEE_29")

    return parser.parse_args()


def load_motif(file_meme_xml, motif_index):
    meme_record = motifs.parse(open(file_meme_xml), 'meme')
    motif = meme_record[motif_index]
    motif.pseudocounts = 1
    return motif


def calculate_background(file):
    sequences = SeqIO.parse(file, "fasta")
    all_sequences = "".join(str(record.seq) for record in sequences)
    nucleotide_frequency = Counter(all_sequences)
    total_nucleotides = sum(nucleotide_frequency.values())
    for nucleotide, count in nucleotide_frequency.items():
        nucleotide_frequency[nucleotide] = count / total_nucleotides
    return nucleotide_frequency


def sm_score(motif, seq):
    q = motif.background
    pwm = np.array([motif.pwm[base] for base in "ACGT"])

    def normalization_factor():
        n = motif.num_occurrences
        return (n + 1) / (n + 4) * log(n + 1) - log(n + 4) - 1 / (n + 4) * sum(log(q[b]) for b in "ACGT") - n / (
                    n + 4) * log(min(q.values()))

    a = normalization_factor()

    def relative_entropy(i):
        return np.sum(pwm[:, i] * np.log(pwm[:, i] / np.array([q[b] for b in "ACGT"]))) / a

    max_score = -float("inf")
    max_index = 0
    motif_len = motif.length
    seq_len = len(seq)

    for i in range(seq_len - motif_len + 1):
        sub_seq = seq[i:i + motif_len]
        sub_seq_indices = np.array([np.where(np.array(list("ACGT")) == base)[0][0] for base in sub_seq])
        score = np.sum(
            [relative_entropy(j) * log(pwm[sub_seq_indices[j], j] / q[sub_seq[j]]) for j in range(motif_len)])
        if score > max_score:
            max_score = score
            max_index = i

    return max_score, seq[max_index:max_index + motif_len], max_index - seq_len


def hamming_distance(s1, s2):
    return np.sum(np.array(list(s1)) != np.array(list(s2)))


def ortholog_sequences(pid, intergen_dir, genome_ref):
    seqs = []
    for file in intergen_dir.iterdir():
        if genome_ref not in file.name:
            for record in SeqIO.parse(file, "fasta"):
                if pid in record.description:
                    seqs.append(record.seq)
    return seqs


def refined_score(motif, pid, intergenic_seq, intergen_dir, genome_ref):
    orthologs = ortholog_sequences(pid, intergen_dir, genome_ref)
    sm_result = sm_score(motif, intergenic_seq)
    s = str(sm_result[1])
    motif_len = motif.length

    ortholog_scores = [
        (motif_len - hamming_distance(s, sm_score(motif, seq)[1])) / motif_len * sm_score(motif, seq)[0]
        for seq in orthologs if len(seq) >= motif_len and "N" not in seq
    ]
    avg_ortholog_score = np.mean(ortholog_scores) if ortholog_scores else 0
    return sm_result[0], sm_result[0] + avg_ortholog_score, s, sm_result[2]


def extract_coding_seq(genome_file):
    for record in SeqIO.parse(genome_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                yield feature.location.extract(record.seq), feature.qualifiers["locus_tag"][0]


def main():
    args = parse_arguments()

    meme_xml = args.meme_xml
    motif_index = args.motif_index
    output_file = args.output_file
    intergen_dir = Path(args.intergen_dir)
    genome_ref = args.genome_ref

    command_line = ET.parse(meme_xml).find("model").find("command_line").text
    motif = load_motif(meme_xml, motif_index)

    for intergen_file in intergen_dir.iterdir():
        if genome_ref in intergen_file.name:
            intergen_path = str(intergen_file)
            break

    motif.background = calculate_background(intergen_path)

    data = [command_line]
    start_time = time.time()
    for record in SeqIO.parse(intergen_path, "fasta"):
        if len(record.seq) < len(motif) or "N" in record.seq:
            continue
        start = record.description.find("WP_")
        end = record.description.find("'", start)
        pid = record.description[start:end]
        result = refined_score(motif, pid, record.seq, intergen_dir, genome_ref)
        data.append(f"{record.id}\t{result[0]}\t{result[1]}\t{result[2]}\t{result[3]}")


    with open(output_file, "w") as f:
        f.writelines("\n".join(data))
    end_time = time.time()
    print(f"Execution time: {end_time - start_time} seconds")


if __name__ == "__main__":
    main()
