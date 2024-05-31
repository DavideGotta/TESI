
def countseqs(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))


def lor(score):
    num = sum(1 for record in SeqIO.parse(intergen, "fasta") if sm(motivo, record.seq)[0] > score) / countseqs(intergen)
    coding = extract_coding_seq(
        "/home/davide/Desktop/genomiChro/annotati_Refseq/Chroococcidiopsis_sp._CCMEE_29_("
        "cyanobacteria)_GCF_023558375.gbff")
    den = sum(1 for seq in coding if sm(motivo, seq[0]) > score)
    if den == 0:
        return 0
    else:
        return log(num / den)

def extract_coding_seq(genoma: str):
    """
    Estrae tutte le sequenze codificanti dal genoma genbank
    """
    for record in SeqIO.parse(genoma, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                yield feature.location.extract(record.seq), feature.qualifiers["locus_tag"][0]

