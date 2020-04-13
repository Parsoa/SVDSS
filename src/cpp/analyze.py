import json

d = []
with open('master_read.txt') as master_seq:
    line = master_seq.readline()
    tokens = line.split()
    for token in tokens:
        i = int(token)
        if i != 65 and i != 67 and i != 71 and i != 84:
            print i
            d.append(token)
with open('master_read.json', 'w') as json_file:
    json.dump(d, json_file, indent = 4)

