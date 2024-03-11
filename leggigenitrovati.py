import re

def extract_words(filename):
    with open(filename, 'r') as file:
        content = file.read()
    matches = re.findall('gene=(\S+)', content)
    return matches

print('\n'.join(extract_words('/home/davide/Downloads/geni_trovati.txt')))