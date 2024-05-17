file="/home/davide/PycharmProjects/TESI2/intergeniche_RefSeq/Chroococcidiopsis_sp._CCMEE_29_GCF_023558375_intergen.fasta"
from Bio import SeqIO
from collections import Counter

def calculate_nucleotide_frequency(file):
    sequences = SeqIO.parse(file, "fasta")
    all_sequences = "".join(str(record.seq) for record in sequences)
    nucleotide_frequency = Counter(all_sequences)
    total_nucleotides =  sum(nucleotide_frequency.values())
    for nucleotide, count in nucleotide_frequency.items():
        nucleotide_frequency[nucleotide] = count / total_nucleotides
    # Calculate the frequency of each dinucleotide (order 1)
    dinucleotide_frequency = Counter(all_sequences[i:i+2] for i in range(len(all_sequences) - 1))
    # Calculate the total number of dinucleotides
    total_dinucleotides = sum(dinucleotide_frequency.values())
    for dinucleotide, count in dinucleotide_frequency.items():
        dinucleotide_frequency[dinucleotide] = count / total_dinucleotides
    return nucleotide_frequency, dinucleotide_frequency

file = "/home/davide/PycharmProjects/TESI2/intergeniche_RefSeq/Chroococcidiopsis_sp._CCMEE_29_GCF_023558375_intergen.fasta"
output_file = "/home/davide/Desktop/Chroococcidiopsis_sp._CCMEE_29_GCF_023558375_intergen_background_model.txt"
with open(output_file, "w") as output:
    nucleotide_frequency, dinucleotide_frequency = calculate_nucleotide_frequency(file)
    print("#    order 0")
    output.write("#    order 0\n")
    for nucleotide, frequency in nucleotide_frequency.items():
        output.write(f"{nucleotide}\t{frequency:.3e}\n")
        print(f"{nucleotide}\t{frequency:.3e}")
    print("#    order 1")
    output.write("#    order 1\n")
    for dinucleotide, frequency in dinucleotide_frequency.items():
        print(f"{dinucleotide}\t{frequency:.3e}")
        output.write(f"{dinucleotide}\t{frequency:.3e}\n")