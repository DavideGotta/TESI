with open("/home/davide/Downloads/KEGGCCMEE.html") as f:
    lines = f.readlines()

#use beautifulsoup to parse the html
from bs4 import BeautifulSoup
soup = BeautifulSoup(open("/home/davide/Downloads/KEGGCCMEE.html"), 'html.parser')
#find all a tags and all b tags
text=soup.get_text()
text=text[:text.find("Organismal Systems")]
lines = text.split("\n")
diz = {}
for i,line in enumerate(lines):
    if line.startswith("0"):
        if "(" in line:
            current_key = line[:line.index("(")].strip()
        else:
            current_key = line
            #add to key lines[i+j] until lines[i+j] starts has (
            j = 1
            while not " (" in lines[i+j]:
                current_key += lines[i+j]
                j += 1
        diz[current_key] = []

    elif line.startswith("WP_"):
        proteins = line.split(", ") if ", " in line else [line]
        diz[current_key].extend(proteins)

import pandas as pd

# Get a list of all unique protein IDs
all_proteins = set(protein for proteins in diz.values() for protein in proteins)

# Initialize an empty list to hold the dataframes
dfs = []

# Iterate over each protein in the set of all proteins
for protein in all_proteins:
    # Create a dictionary where the keys are the pathways and the values are 1 if the protein is present in the pathway and 0 otherwise
    data = {pathway: 1 if protein in proteins else 0 for pathway, proteins in diz.items()}
    # Convert the dictionary into a DataFrame and append it to the list
    df = pd.DataFrame(data, index=[protein])
    dfs.append(df)

# Concatenate all the dataframes into a single dataframe
result = pd.concat(dfs)
