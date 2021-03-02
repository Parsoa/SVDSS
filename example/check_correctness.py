import sys

def load_solution(path):
    n = 0
    s_entries = {}
    with open(path) as s_file:
        line = s_file.readline()
        while line:
            if n % 4 == 0:
                count = line.split('#')[1]
            if n % 4 == 1:
                s_entries[line.strip()] = count
            n += 1
            line = s_file.readline()
    return s_entries

a = load_solution(sys.argv[1])
b = load_solution(sys.argv[2])

for key in a:
    if key in b and a[key] == b[key]:
        pass
    else:
        print('Input mismatch.')
        exit()

print('Input match.')
