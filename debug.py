file="/home/davide/Downloads/genomiChro/blast_results_all/best_hits/regLexA_VS_Chroococcidiopsis_thermalis_PCC_7203_(cyanobacteria)_GCF_000317125.1_best_hits.txt"
with open(file, 'r') as f:
    content = f.read()
count=0
for i,line in enumerate(content.split("\n")):
    if line.startswith(">"):
        count+=1
print(count)