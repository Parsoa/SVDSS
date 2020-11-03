#!/usr/bin/python

import os
import sys

n = 0
m = 0

out_file = open(sys.argv[3], 'w')

diff_entries = []

d = 0
with open(sys.argv[2]) as diff_file:
    diff_line = diff_file.readline()
    diff_line = diff_file.readline()
    while diff_line:
        diff_entries.append(diff_line.strip())
        diff_line = diff_file.readline()

i = 0
print('Loaded', len(diff_entries), ' diff entries.')

with open(sys.argv[1]) as fasta_file:
    fasta_header = fasta_file.readline()
    fasta_sequence = fasta_file.readline()
    while fasta_header:
        prev = 0
        read = ''
        offset = 0
        read_id = fasta_header.strip()[1:]
        while i < len(diff_entries) and diff_entries[i].split()[0] == read_id:
            tokens = diff_entries[i].split()
            position = int(tokens[1]) - 1
            target = tokens[3]
            if '-' in target:
                read += fasta_sequence[prev: position + offset]
                read += tokens[2]
                prev = position + offset
                offset -= 1
            elif '+' in target:
                read += fasta_sequence[prev: position + offset]
                prev = position + offset + 1
                offset += 1
            else:
                read += fasta_sequence[prev: position + offset]
                read += tokens[2]
                prev = position + offset + 1
                offset = offset
            i += 1
        read += fasta_sequence[prev:]
        out_file.write('@' + read_id + '\n')
        out_file.write(read)
        out_file.write('+' + read_id + '\n')
        out_file.write('~' * (len(read) - 1) + '\n')
        m += 1
        fasta_header = fasta_file.readline()
        fasta_sequence = fasta_file.readline()
        n += 1
        if n % 1000 == 0:
            print('In:', n, 'reads.')
            print('Out:', m, 'reads.')
