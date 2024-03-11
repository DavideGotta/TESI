from Bio.motifs import meme

with open("/home/davide/Downloads/meme.xml") as f:

    record = meme.read(f)

for motif in record:
    for instance in motif.instances:
        print(instance.motif_name, instance.sequence_name, instance.strand, instance.pvalue)