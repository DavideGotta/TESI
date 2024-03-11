import re
geni = []
with open("/home/davide/Downloads/geni_trovati(2).txt", "r") as f:
    #find lines duplicates
    lines = f.readlines()
    for line in lines:
        gene = re.search(r'gene=(\w+)', line).group(1)
        geni.append(gene)
from collections import Counter
for gene in sorted(Counter(geni), key=Counter(geni).get, reverse=True):
    print(f"{gene} {Counter(geni)[gene]}")