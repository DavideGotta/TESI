# Purpose: Parse KEGG BRITE html file to extract the hierarchy of the brites
#use beautifulsoup to parse the html
from bs4 import BeautifulSoup
import pandas as pd
def parsekeggbrite(file):
    soup = BeautifulSoup(open(file), 'html.parser')
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
    all_proteins = set(protein for proteins in diz.values() for protein in proteins)
    dfs = []
    for protein in all_proteins:
        data = {pathway: 1 if protein in proteins else 0 for pathway, proteins in diz.items()}
        df = pd.DataFrame(data, index=[protein])
        dfs.append(df)
    result = pd.concat(dfs)
    return result

def parsekeggpathway(file):
    soup = BeautifulSoup(open(file), 'html.parser')
    text = soup.get_text()
    text = text[:text.find("Organismal Systems")]
    lines = text.split("\n")
    diz = {}
    for i, line in enumerate(lines):
        if line.startswith("0"):
            if "(" in line:
                current_key = line[:line.index("(")].strip()
            else:
                current_key = line
                j = 1
                while not " (" in lines[i + j]:
                    current_key += lines[i + j]
                    j += 1
            diz[current_key] = []

        elif line.startswith("WP_"):
            if current_key not in diz:
                continue
            proteins = line.split(", ") if ", " in line else [line]
            diz[current_key].extend(proteins)

    import pandas as pd

    # Get a list of all unique protein IDs
    all_proteins = set(protein for proteins in diz.values() for protein in proteins)
    # Initialize an empty list to hold the dataframes
    dfs = []
    # Iterate over each protein in the set of all proteins
    for protein in all_proteins:
        # Create a dictionary where the keys are the brites and the values are 1 if the protein is present in the pathway and 0 otherwise
        data = {pathway: 1 if protein in proteins else 0 for pathway, proteins in diz.items()}
        # Convert the dictionary into a DataFrame and append it to the list
        df = pd.DataFrame(data, index=[protein])
        dfs.append(df)

    # Concatenate all the dataframes into a single dataframe
    result = pd.concat(dfs)
    return result
from pprint import pprint
file=("/home/davide/Downloads/KEGGBrite.html")
soup=BeautifulSoup(open(file), 'html.parser')
li_tags = soup.find_all("li")
for i, li in enumerate(li_tags):
    if "Enzymes" in str(li):
        break
li_tags = li_tags[i:]

brites = {}
#now iterate over the li tags and find the dd tags inside them
for li in li_tags:
    dd_tags = li.find_all("dd")
    brite=str(li.contents[1])[:-1].strip()
    diz = {}
    for dd in dd_tags:
        contents = dd.contents
        contents=[c for c in contents if c.name != 'br']
        contents=[str(content).replace('\xa0', '') for content in contents]
        contents=[str(c) for c in contents]
        for i, e in enumerate(contents):
            e=e.split(":")
            pids=e[1].split(", ")

            key=e[0]
            if key not in diz:
                diz[key]=pids
            elif pids!=['']:
                diz[key].extend(pids)
    brites[brite]=diz

for brite in brites:
    diz=brites[brite]
    for key in diz:
        diz[key]=len(set(diz[key]))

pprint(brites)
#delete first elements unitl one has Metabolic patwhays in text
