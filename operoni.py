
# Specify the input file
file="/home/davide/Desktop/genomiChro/operons/list_of_operons_3500660"
#parse the file
with open(file, 'r') as file:
    content = file.read()
    #split if there is number as first character in a line
    lines = content.strip().split("\n")
    #strip the lines
    #if the line starts with a number, that line and all the following lines until the next number are part of the same operon

    operons = {}
    for line in lines[1:]:
        if line[0].isdigit():
            header = line.split(" ")[0]
            print(header)
            operons[header] = []
            if header=='3807':
                break
        else:
            operons[header].append(line.strip().split("\t"))


    #print the number of operons

operons_filtered = {k: v for k, v in operons.items() if len(v) > 1}

#write the operons filtered to a new file
with open("operoni", 'w') as file:
    for k,v in operons_filtered.items():
        inter=v[0][3]+" "+v[-1][4]+" "+v[0][5]
        print(k,inter)
        file.write(inter+"\n")
        #file.write("\t".join(line) + "\n")