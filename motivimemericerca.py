import os
from Bio import SeqIO
file="/home/davide/Downloads/fimo.tsv"
file2="/home/davide/Downloads/genomiChro/gffs/Chroococcidiopsis_sp._CCMEE_29_cyanobacteria_GCF_023558375.1.gff"
#open the gff file and read the lines
i=0
geni= []
with open(file, "r") as f:
    lines = f.readlines()
    #for each line, split the line by tab
    for line in lines:
        line = line.split("\t")
        if len(line)<=2:    #if the line is empty, skip
            continue

        with open(file2, "r") as f2:
            lines2 = f2.readlines()
            for line2 in lines2:
                if line[2] in line2 and "gene=" in line2:
                    geni.append(line2.split("\t")[-1].split(";")[-4].split("=")[-1])
                    i+=1
                    break
print(i)
print(geni)
dir="/home/davide/Downloads/fimo/"
gff_dir="/home/davide/Downloads/genomiChro/gffs/"

import os
#open all gff files in the directory dir
for file in os.listdir(dir):
    with open(os.path.join(dir,file), "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.split("\t")
            if len(line)<=2:    #if the line is empty, skip
                continue
            for file2 in os.listdir(gff_dir):
                if file[4:-12] in file2:
                    with open(os.path.join(gff_dir,file2), "r") as f2:
                        lines2 = f2.readlines()
                        for line2 in lines2:
                            if line[2] in line2 and "gene=" in line2:
                                geni.append(line2.split("\t")[-1].split(";")[-4].split("=")[-1])
                                break
from collections import Counter
with open("genifimo.txt", "w") as f:
    for gene in sorted(Counter(geni), key=Counter(geni).get, reverse=True):
        f.write(f"{gene} {Counter(geni)[gene]}\n")