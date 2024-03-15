import re  # Import the re module for regular expressions
from Bio import SeqIO  # Import the SeqIO module from the Biopython package
from Bio.SeqRecord import SeqRecord  # Import the SeqRecord module from the Biopython package
from Bio.Seq import Seq  # Import the Seq module from the Biopython package
from Bio.SeqFeature import SeqFeature, FeatureLocation  # Import the SeqFeature module from the Biopython package
from Bio.Seq import MutableSeq  # Import the MutableSeq module from the Biopython package
def geneDentroOperone(gene_start, fil, ranges):
    """
    :param gene_start:  start of the gene 
    :param fil:     strand of the gene
    :param ranges:  list of tuples with start and end of the genes in the operon
    :return:    True if the gene is inside an operon, False otherwise
    """
    for line in ranges:
        start, end, strand = line.split()
        start, end = int(start), int(end)
        if strand == fil:
            if start < gene_start < end:
                return True
    return False


def geneInizioOperone(gene_start, fil, ranges):
    """
    :param gene_start:  start of the gene
    :param fil:     strand of the gene
    :param ranges:  list of tuples with start and end of the genes in the operon
    :return:    True if the gene is at the start of an operon, False otherwise
    """
    for line in ranges:
        start, end, strand = line.split()
        start, end = int(start), int(end)
        if strand == fil and start == gene_start:
            return True
    return False


# Define a function to find intergenic sequences
def sequenze_intergeniche(genoma, geni, output):
    geni, e_score, perc_id, queries = zip(*geni)  # Unzip the list of genes of interest
    intergeniche = []  # List to store intergenic sequences
    # Read the GenBank file
    genomi = SeqIO.parse(genoma, "genbank")
    # Loop through the features in the GenBank file
    for genoma in genomi:
        cds_list_plus = []  # List to store gene info in + direction
        cds_list_minus = []  # List to store gene info in - direction
        for feature in genoma.features:
            # Check if the feature is a CDS
            if feature.type == "CDS":
                mystart = feature.location.start  # Get the start of the CDS
                myend = feature.location.end  # Get the end of the CDS
                # Check the strand of the CDS and add the info to the appropriate list
                if feature.strand == -1:
                    cds_list_minus.append((mystart, myend, feature.qualifiers["locus_tag"][0],
                                           feature.qualifiers["product"][0], feature.qualifiers.get("gene", [' '])[0]))
                else:
                    cds_list_plus.append((mystart, myend, feature.qualifiers["locus_tag"][0],
                                          feature.qualifiers["product"][0], feature.qualifiers.get("gene", [' '])[0]))
        # Loop through the genes in + direction
        for i, n in enumerate(cds_list_plus):
            last_end = cds_list_plus[i - 1][1]  # Get the end of the last gene
            this_start = n[0]  # Get the start of the current gene
            # Check if the start of the current gene is in the list of starts
            if cds_list_plus[i][2] in geni:
                indice = geni.index(cds_list_plus[i][2])
                # Get the intergenic sequence
                if last_end > this_start:
                    last_end = 0
                intergene_seq = genoma.seq[last_end:this_start]

                # If the intergenic sequence is longer than 300, only take the last 300 bases
                if len(intergene_seq) >= 300 or len(intergene_seq) == 0:
                    intergene_seq = genoma.seq[this_start - 300:this_start]
                    last_end = this_start - 300

                # Add the intergenic sequence to the list as a SeqRecord
                intergeniche.append(
                    SeqRecord(
                        intergene_seq,
                        id=n[2],
                        description=f"product = {cds_list_plus[i][3]}, gene = {cds_list_plus[i][4]} similar to {queries[indice]} with E-value {e_score[indice]} and identity {perc_id[indice]}% {last_end + 1}-{this_start} +"
                    )
                )
        # Loop through the genes in - direction
        for i, n in enumerate(cds_list_minus):
            if i == len(cds_list_minus) - 1:
                i = -1
            next_start = cds_list_minus[i + 1][0]  # Get the start of the next gene
            this_end = n[1]  # Get the end of the current gene
            # Check if the start of the current gene is in the list of starts
            if cds_list_minus[i][2] in geni:
                indice = geni.index(cds_list_minus[i][2])
                # Get the intergenic sequence and reverse complement it
                if this_end > next_start:
                    next_start = len(genoma.seq)

                intergene_seq = genoma.seq[this_end:next_start]

                # If the intergenic sequence is longer than 300, only take the last 300 bases
                if len(intergene_seq) >= 300 or len(intergene_seq) == 0:
                    intergene_seq = genoma.seq[this_end:this_end + 300]
                    next_start = this_end + 300
                intergene_seq = intergene_seq.reverse_complement()
                # Add the intergenic sequence to the list as a SeqRecord

                intergeniche.append(
                    SeqRecord(
                        intergene_seq,
                        id=n[2],
                        description=f"product = {cds_list_minus[i][3]}, gene = {cds_list_minus[i][4]} similar to {queries[indice]} with E-value {e_score[indice]} and identity {perc_id[indice]}% {this_end}-{next_start + 1} -"
                    )
                )
    SeqIO.write(intergeniche, output, "fasta")


def sequenze_intergenichetutte(input, output, max_len=300,op=False):
    intergeniche = []  # Lista per memorizzare le sequenze intergeniche
    genoma = SeqIO.parse(input, "genbank")

    # Loop attraverso le features nel file GenBank
    for sequenze in genoma:
        # sequenze sta per sequenzefold perchè il file genbank potrebbe non contenere un solo cromosoma circolare
        cds_list_plus = []  # Lista per memorizzare le informazioni sui geni in direzione +
        cds_list_minus = []  # Lista per memorizzare le informazioni sui geni in direzione -
        for feature in sequenze.features:
            # Controlla se la feature è un CDS
            if feature.type == "CDS":
                mystart = feature.location.start  # Ottieni l'inizio del CDS
                myend = feature.location.end  # Ottieni la fine del CDS
                mystrand = '-' if feature.location.strand == -1 else '+'
                mystart, myend = int(mystart), int(myend)
                # Check the strand of the CDS and add the info to the appropriate list
                if op:
                    with open("operoniCCMEE", "r") as f:
                        content = f.readlines()
                    if not geneDentroOperone(mystart + 1, mystrand, content):
                        if mystrand == '-':
                            cds_list_minus.append((mystart, myend, feature.qualifiers["locus_tag"][0],
                                                   feature.qualifiers["product"][0],
                                                   feature.qualifiers.get("gene", [' '])[0],
                                                   geneInizioOperone(mystart + 1, mystrand, content)))
                        else:
                            cds_list_plus.append((mystart, myend, feature.qualifiers["locus_tag"][0],
                                                  feature.qualifiers["product"][0],
                                                  feature.qualifiers.get("gene", [' '])[0],
                                                  geneInizioOperone(mystart + 1, mystrand, content)))
                else:
                    if mystrand == '-':
                        cds_list_minus.append((mystart, myend, feature.qualifiers["locus_tag"][0],
                                               feature.qualifiers["product"][0],
                                               feature.qualifiers.get("gene", [' '])[0]))
                    else:
                        cds_list_plus.append((mystart, myend, feature.qualifiers["locus_tag"][0],
                                              feature.qualifiers["product"][0],
                                              feature.qualifiers.get("gene", [' '])[0]))

        # Loop attraverso i geni in direzione +
        for i, n in enumerate(cds_list_plus):

            last_end = 0 if i == 0 else cds_list_plus[i - 1][1]
            this_start = n[0]  # Ottieni l'inizio del gene corrente
            intergene_seq = sequenze.seq[last_end:this_start]  # Ottieni la sequenza intergenica
            if len(intergene_seq) >= max_len:  # Se la sequenza intergenica è più lunga di max_len, prendi solo le ultime max_len basi(default 300)
                intergene_seq = sequenze.seq[this_start - max_len:this_start]
                last_end = this_start - max_len  # Aggiorna quello che sarà l'inizio della sequenza intergenica
            # Aggiungi la sequenza intergenica alla lista come SeqRecord
            description = f"product = {cds_list_plus[i][3]}, gene = {cds_list_plus[i][4]} {last_end + 1}-{this_start} +"
            if op and n[-1]:
                description = "OPERON " + description
            if len(intergene_seq) == 0:
                continue
            intergeniche.append(
                SeqRecord(
                    intergene_seq,
                    id=n[2],  # id del gene(locus_tag)
                    description=description
                )
            )
        # Loop attraverso i geni in direzione -
        for i, n in enumerate(cds_list_minus):
            next_start = len(sequenze.seq) if i + 1 == len(cds_list_minus) else cds_list_minus[i + 1][
                0]  # Ottieni l'inizio del prossimo gene(gene precedente per la direzione -)
            this_end = n[1]  # Ottieni la fine del gene corrente(inizio del gene per la direzione -)
            intergene_seq = sequenze.seq[this_end:next_start]  # Ottieni la sequenza intergenica
            if len(intergene_seq) >= max_len:
                # Se la sequenza intergenica è più lunga di max_len, prendi solo le prime max_len basi(default 300)
                intergene_seq = sequenze.seq[this_end:this_end + max_len] if this_end + max_len < len(
                    sequenze.seq) else sequenze.seq[this_end:]
                next_start = this_end + max_len
            intergene_seq = intergene_seq.reverse_complement()  # Complemento inverso della sequenza intergenica
            description = f"product = {cds_list_minus[i][3]}, gene = {cds_list_minus[i][4]} {this_end}-{next_start + 1} -"
            if op and n[-1]:
                description = "OPERON " + description
            # Aggiungi la sequenza intergenica alla lista come SeqRecord
            if len(intergene_seq) == 0:
                continue
            intergeniche.append(
                SeqRecord(
                    intergene_seq,
                    id=n[2],  # id del gene(locus_tag)
                    description=description
                )
            )

    # Scrivi le sequenze intergeniche nel file di output
    SeqIO.write(intergeniche, output, "fasta")

