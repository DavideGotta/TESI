# Purpose: Parse KEGG BRITE html file to extract the hierarchy of the pathways
#use beautifulsoup to parse the html
from bs4 import BeautifulSoup
soup = BeautifulSoup(open("/home/davide/Downloads/KEGGbriteCCMEE29.html"), 'html.parser')
text=soup.get_text()
text=text[text.find("ko01001"):]
lines = text.split("\n")
diz={}

for i,line in enumerate(lines):
    if line.startswith("ko"):
        current_key = line[:line.index("(")].strip()
        diz[current_key] = []
    elif line.startswith("WP"):
        proteins = line.split(",") if "," in line else [line]
        diz[current_key].extend(proteins)
from pprint import pprint

import pandas as pd
all_proteins = set(protein for proteins in diz.values() for protein in proteins)
print(len(all_proteins))
print(all_proteins)
dfs = []
for protein in all_proteins:
    data = {pathway: 1 if protein in proteins else 0 for pathway, proteins in diz.items()}
    df = pd.DataFrame(data, index=[protein])
    dfs.append(df)
result = pd.concat(dfs)

print(result)
sum=0
#delete lines empty
lines = [line for line in lines if line]
for line in lines[:-5]:
    num=int(line[line.find("(")+1:line.find(")")]) if "(" in line else 0
    sum+=num


print(sum)
print(lines)
for key, value in diz.items():
    print(key, len(value))