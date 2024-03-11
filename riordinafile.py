import re
import collections
file="/home/davide/Downloads/genomiChro/intergeniche_best_hits/motivitrovati/joined.txt"
# Read the filefile
with open(file, 'r') as f:
    content = f.read()
sep="*"*100
# Split the content into blocks
blocks = content.split(sep)
strings = []
lines = content.split("\n")
# Iterate over each line in the file
for i in range(len(lines)):
    # If the line contains a string inside [] square brackets
    match = re.search(r'\'(.*?)\'', lines[i])
    if match and i + 2 < len(lines) and lines[i + 2].startswith('Trovato'):

        # Add the string to the list
        strings.append(match.group(1))

# Count how many times each unique string occurs
counts = collections.Counter(strings)

#sort for biggest number of occurences
sorted_counts = sorted(counts.items(), key=lambda x: x[1])
diz = dict(sorted_counts)
diz['']=0
# Print the counts
print(sorted_counts)
# Define a function to find the string inside square brackets
def find_string_in_brackets(block):
    match = re.search(r'\'(.*?)\'', block)
    if match:
        return match.group(1)
    else:
        return ''
from collections import Counter

# Count the frequency of each string inside square brackets
counter = Counter(find_string_in_brackets(block) for block in blocks)

# Sort the blocks based on the frequency and the string inside square brackets
sorted_blocks = sorted(blocks, key=lambda block: (-diz[find_string_in_brackets(block)], find_string_in_brackets(block)))

# Join the sorted blocks
sorted_content = sep.join(sorted_blocks)

# Write the sorted content back to the file
with open(file[:-4]+"sorted.txt", 'w') as f:
    f.write(sorted_content)
