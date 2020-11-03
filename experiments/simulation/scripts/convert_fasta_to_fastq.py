#!/usr/bin/python

with open('reads_edited.fa') as fasta_file:
    with open('reads.fastq', 'w') as fastq_file:
        header = fasta_file.readline()
        seq = fasta_file.readline()
        while header:
            fastq_file.write('@' + header[1:])
            fastq_file.write(seq)
            fastq_file.write('+' + header[1:])
            fastq_file.write('I' * (len(seq) - 1) + '\n')
            header = fasta_file.readline()
            seq = fasta_file.readline()

