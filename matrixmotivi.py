import pandas as pd
file="/home/davide/Desktop/genomiChro/intergeniche_tutte/motivounicobest/joined.txt"
with open(file, "r") as f:
    #if line starts with Risultati per il gene: extract the gene name
    lines = f.readlines()

#create dataframe with columns 'HLAHILNB', 'LHCOABFD', 'JHIDLJMF', 'AIHNMJLC', 'GOENLOBA', 'FNLHCIEN', 'GBPHKEKP', 'GOHPBNEF', 'MGCKOBKL'
df = pd.DataFrame(columns=['HLAHILNB', 'LHCOABFD', 'JHIDLJMF', 'AIHNMJLC', 'GOENLOBA', 'FNLHCIEN', 'GBPHKEKP', 'GOHPBNEF', 'MGCKOBKL'])
#iterate over the genes and add a row for each gene if the gene is not already in the dataframe
i=0
motivi=0
for line in lines:
    if i==10:
        break
    if line.startswith("Risultati per il gene:"):
        i+=1
        if motivi!=0:
            df.loc[id,name_col]=motivi

        id = line[line.find('[')+2:line.find(']')-1].split(" ")[0]
        name_col=line.split(" ")[4].split("_")[0]
        motivi=[]
    else:
        if "il motivo" in line:
            motivi.append(line.split("-")[1].split(" ")[2]+"("+line.split("-")[2][:-1]+")")

print(df)