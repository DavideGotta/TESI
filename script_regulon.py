import math  # Importa il modulo math per fornire funzioni matematiche
import random  # Importa il modulo random per generare numeri casuali
from Bio import SeqIO  # Importa la classe SeqIO dal modulo Biopython per la gestione delle sequenze
from Bio.Seq import Seq  # Importa la classe Seq dal modulo Biopython per la manipolazione delle sequenze
from Bio.SeqRecord import SeqRecord  # Importa la classe SeqRecord dal modulo Biopython per rappresentare le sequenze


# da aggiungere una funzione che calcoli la distanza di Hamming e darla come ortholog_distances a calculate_refined_score
# ortholog_distances è una lista di distanze di Hamming tra la sequenza data (che è la sequenza per la quale stiamo calcolando il punteggio raffinato) e le sequenze ortologhe.

def calculate_q_bases(inter_TU_regions):
    # Funzione per calcolare le frequenze delle basi nei regioni inter-TU
    q_bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0}  # Inizializza un dizionario per le frequenze delle basi
    total_bases = 0  # Inizializza il conteggio totale delle basi

    # Calcola le frequenze delle basi per tutte le regioni inter-TU
    for region in inter_TU_regions:
        for base in region:
            if base in q_bases:
                q_bases[base] += 1  # Incrementa il conteggio della base
                total_bases += 1  # Incrementa il conteggio totale delle basi

    # Calcola le frequenze relative dividendo per il numero totale di basi
    for base in q_bases:
        q_bases[base] /= total_bases  # Calcola la frequenza relativa

    return q_bases  # Restituisce le frequenze delle basi


def SM_score(sequence, profile, q_bases, n):
    # Funzione per calcolare il punteggio SM per una sequenza
    max_score = float('-inf')  # Inizializza il punteggio massimo come infinito negativo
    profile_length = len(profile)  # Lunghezza del profilo
    sequence_length = len(sequence)  # Lunghezza della sequenza

    # Calcola il fattore di normalizzazione
    a = normalization_factor(n, q_bases)

    # Calcola il punteggio SM per ogni sottostringa h della sequenza
    for i in range(sequence_length - profile_length + 1):
        substring = sequence[i:i + profile_length]  # Estrae la sottostringa della sequenza
        score = 0  # Inizializza il punteggio per la sottostringa

        # Calcola il contributo logaritmico per ciascuna posizione nel profilo
        for j in range(profile_length):
            base = substring[j]  # Base nella posizione j della sottostringa
            score += math.log(profile[j][base] / q_bases[base])  # Aggiunge il contributo al punteggio

        # Sottrae il fattore di normalizzazione
        score -= a
        # Aggiorna il punteggio massimo se necessario
        if score > max_score:
            max_score = score  # Aggiorna il punteggio massimo se necessario

    return max_score  # Restituisce il punteggio massimo


def information_i(profile_column, q_bases, n):
    # Funzione per calcolare il contenuto informativo I_i (equazione 2) di una colonna del profilo
    info_content = 0  # Inizializza il contenuto informativo a zero
    a = normalization_factor(n, q_bases)  # Calcola il fattore di normalizzazione

    # Calcola il contenuto informativo della colonna
    for base, frequency in profile_column.items():
        info_content += frequency * math.log(frequency / q_bases[base])  # Aggiunge il contributo informativo

    # Normalizza dividendo per il fattore di normalizzazione
    info_content /= a

    return info_content  # Restituisce il contenuto informativo

    # La funzione information_i accetta tre argomenti: profile_column, q_bases, e n. Questi sono:
    # profile_column: un dizionario che rappresenta una colonna del profilo. Le chiavi sono le basi {A, C, G, T}, e i valori sono le frequenze di quelle basi nella colonna.
    # q_bases: un dizionario che rappresenta le frequenze delle basi nelle regioni inter-TU aggregate. Le chiavi sono le basi {A, C, G, T}, e i valori sono le frequenze di quelle basi.
    # n: il numero di siti di legame per costruire il profilo M.
    # Calcolo del contenuto informativo:
    # La funzione calcola il contenuto informativo della colonna del profilo utilizzando la formula fornita:
    # Per ogni base b nelle basi {A, C, G, T}, la funzione itera attraverso le basi presenti nella colonna del profilo e calcola p(i,b).  Questo valore rappresenta il contributo informativo della base b nella colonna del profilo.
    # Questi contributi informativi vengono sommati per tutte le basi e divisi per il fattore di normalizzazione a.
    # Il risultato finale è il contenuto informativo della colonna del profilo.
    # Restituzione del risultato:
    # La funzione restituisce il contenuto informativo calcolato per la colonna del profilo.


def normalization_factor(n, q_bases):
    # Funzione per calcolare il fattore di normalizzazione
    a = (math.log(n + 1) - math.log(n + 4)) / (n + 4 - n / (n + 4))  # Calcola il fattore di normalizzazione a
    min_background_frequency = min(
        q_bases.values())  # Trova la frequenza minima tra le basi nelle regioni inter-TU aggregate
    b_sum = sum([math.log(frequency) - (n / (n + 4)) * math.log(min_background_frequency) for frequency in
                 q_bases.values()])  # Calcola la somma dei logaritmi delle frequenze delle basi nelle regioni inter-TU aggregate
    return b_sum - a  # Restituisce la differenza tra la somma dei logaritmi e il fattore di normalizzazione


def calculate_refined_score(sequence, profile_scores, ortholog_sequences, ortholog_distances):
    # Funzione per calcolare il punteggio raffinato
    base_score = sum(profile_scores)  # Calcola il punteggio base come la somma dei punteggi del profilo
    max_score = float('-inf')  # Inizializza il punteggio massimo come infinito negativo

    # Calcola il punteggio raffinato
    for ortholog_sequence, ortholog_distance in zip(ortholog_sequences, ortholog_distances):
        # Calcola il punteggio parziale per l'ortologo
        partial_score = sum(profile_scores[i] - ortholog_distance[i] for i in
                            range(len(profile_scores)))  # Calcola il punteggio parziale per l'ortologo
        # Aggiorna il punteggio massimo
        max_score = max(max_score, partial_score)  # Aggiorna il punteggio massimo se necessario

    # Aggiorna il punteggio raffinato
    refined_score = base_score + max_score  # Aggiorna il punteggio raffinato

    return refined_score  # Restituisce il punteggio raffinato


def extract_coding_sequences(genome_file, gff_file):
    # Funzione per estrarre le sequenze codificanti dal file del genoma
    coding_regions = {}  # Dizionario per memorizzare le coordinate delle regioni codificanti per ciascun contig

    # Lettura del file GFF per estrarre le coordinate delle regioni codificanti
    with open(gff_file, 'r') as gff:
        for line in gff:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                if parts[2] == 'CDS':  # Se la feature è CDS (Coding Sequence)
                    contig_id = parts[0]  # ID del contig
                    start = int(parts[3])  # Posizione di inizio della CDS
                    end = int(parts[4])  # Posizione di fine della CDS
                    if contig_id not in coding_regions:
                        coding_regions[
                            contig_id] = []  # Crea una nuova lista per le coordinate della CDS per questo contig
                    coding_regions[contig_id].append((start, end))  # Aggiungi le coordinate alla lista

    # Estrazione delle sequenze codificanti dal file del genoma
    print(coding_regions)
    coding_sequences = []
    for record in SeqIO.parse(genome_file, 'fasta'):  # Per ogni record nel file del genoma in formato FASTA
        contig_id = record.id  # ID del contig
        if contig_id in coding_regions:  # Se ci sono regioni codificanti per questo contig
            for start, end in coding_regions[contig_id]:  # Per ogni coppia di coordinate di regioni codificanti
                sequence = record.seq[start - 1:end]  # Estrai la sequenza codificante (1-indexed, quindi sottrai 1)
                print(record)
                coding_sequences.append(SeqRecord(sequence, id=f"{contig_id}_CDS_{start}_{end}",
                                                  description=""))  # Aggiungi la sequenza all'elenco
    # print(coding_sequences)
    return coding_sequences  # Restituisci l'elenco delle sequenze codificanti estratte


def generate_random_sequences(coding_sequences, num_sequences):
    # Funzione per generare sequenze casuali
    random_sequences = []
    for _ in range(num_sequences):  # Per il numero desiderato di sequenze casuali
        random_sequence = ''.join(random.choice('ACGT') for _ in range(len(
            coding_sequences[0].seq)))  # Genera una sequenza casuale della stessa lunghezza delle sequenze codificanti
        random_sequences.append(SeqRecord(Seq(random_sequence), id=f"random_{_ + 1}",
                                          description=""))  # Aggiungi la sequenza casuale all'elenco

    return random_sequences  # Restituisci l'elenco delle sequenze casuali


# Esempio di come utilizzarlo
genome_file = "/home/davide/Desktop/genomiChro/fasta/Chroococcidiopsis_sp_CCMEEE_29_chromosome.fa"  # Percorso del file del genoma
gff_file = "/home/davide/Downloads/sequence.gff3"  # Percorso del file GFF
num_random_sequences = 300  # Numero di sequenze casuali da generare

# Estrai le sequenze codificanti dal genoma basate sul file GFF
coding_sequences = extract_coding_sequences(genome_file, gff_file)
print(coding_sequences)

# Genera sequenze casuali di lunghezza uguale alle sequenze codificanti
random_sequences = generate_random_sequences(coding_sequences, num_random_sequences)


# Ora si possono utilizzare queste sequenze casuali come CU per valutare la significatività statistica delle predizioni.

def calculate_false_positive_rate(IU_lengths, profile, num_random_seqs, threshold):
    # Funzione per calcolare il tasso di falsi positivi
    random_scores = [calculate_score(generate_random_sequence(length), profile) for length in IU_lengths for _ in
                     range(num_random_seqs)]  # Genera punteggi casuali
    num_positive_random = sum(score > threshold for score in random_scores)  # Conta punteggi casuali sopra la soglia
    false_positive_rate = num_positive_random / len(random_scores)  # Calcola il tasso di falsi positivi
    return false_positive_rate  # Restituisce il tasso di falsi positivi


def calculate_log_odds_ratio(IU_scores, CU_scores, threshold):
    # Funzione per calcolare il rapporto di log odds
    p_IU = sum(score > threshold for score in IU_scores) / len(
        IU_scores)  # Frazione di sequenze inter-TU con punteggio > soglia
    p_CU = sum(score > threshold for score in CU_scores) / len(
        CU_scores)  # Frazione di sequenze casuali con punteggio > soglia
    LOR = math.log(p_IU / p_CU)  # Calcola il rapporto di log odds
    return LOR  # Restituisce il rapporto di log odds