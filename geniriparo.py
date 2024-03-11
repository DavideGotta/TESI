file="/home/davide/geniriparo.txt"
import re
#extract all names after a number followed by a point like 1. genC extract genC
i=0
with open(file,"r") as f:
    lines = f.readlines()
    for line in lines:
        try:
            gene = re.search(r'\d+\.\s(\w+)', line).group(1)
            print(gene)
            i+=1
        except:
            pass
print(i)