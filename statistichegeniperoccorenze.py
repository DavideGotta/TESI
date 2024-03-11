file = "/home/davide/Desktop/genomiChro/intergeniche_best_hits/motivitrovati/genilistaperoccorrenze.txt"
import re
from collections import Counter
#extract all numbers after Occorenze = for example Occorenze = 10, count how how many times each number appears
with open(file, "r") as f:
    occ = f.read()
    occ = re.findall(r'Occorrenze = (\d+)', occ)
    occ = [int(i) for i in occ]
    occ = Counter(occ)
    print(occ)
#use matplotlib to plot the results
#use seaborn to make the plot more readable
import matplotlib.pyplot as plt
import numpy as np

plt.bar(occ.keys(), occ.values())
plt.xlabel("Occorrenze")
plt.ylabel("Numero di geni")

# Set x-axis labels from 1 to 10
plt.xticks(np.arange(1, 11))
plt.bar(occ.keys(), occ.values())
plt.xlabel("Numero di genomi in cui compare motivo")
plt.ylabel("Numero di geni")

# Set x-axis labels from 1 to 10
plt.xticks(np.arange(1, 11))

# Save the figure
plt.savefig('output.png')


plt.show()