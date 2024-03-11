def extract_best_hits(input_filename, output_filename):
    with open(input_filename, 'r') as file:
        lines = file.readlines()
    with open(output_filename, 'w') as outfile:
        for i in range(len(lines)):
            if lines[i].startswith('Query='):
                outfile.write(lines[i]+lines[i+1])
            if 'Sequences producing significant alignments' in lines[i]:
                for j in range(i+1, len(lines)):
                    if lines[j].startswith('>'):
                        k = 0
                        while True:  # While the line is not blank
                            outfile.write(lines[j + k])
                            if "Identities" in lines[j + k]:
                                break
                            k += 1
                        outfile.write('\n')
                        break

#write_next_line('/home/davide/Downloads/TS821/genilexAvsTS821', '/home/davide/Downloads/TS821/genisimili')
import re

def extract_words(filename):
    with open(filename, 'r') as file:
        content = file.read()
    matches = re.findall(r'>(\S+)', content)
    return matches

geni=extract_words('/home/davide/Downloads/TS821/genisimili')
