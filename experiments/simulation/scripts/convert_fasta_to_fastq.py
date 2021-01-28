#!/usr/bin/python
import sys

with open(sys.argv[1]) as fasta_file:
    with open(sys.argv[2], 'w') as fastq_file:
        header = fasta_file.readline()
        seq = fasta_file.readline()
        while header:
            fastq_file.write('@' + header[1:])
            fastq_file.write(seq)
            fastq_file.write('+' + header[1:])
            fastq_file.write('I' * (len(seq) - 1) + '\n')
            header = fasta_file.readline()
            seq = fasta_file.readline()

