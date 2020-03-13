from __future__ import print_function

import io
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import random
import traceback

from fulminata import (
    config,
)

from fulminata.logger import *

# ============================================================================================================================ #
# kmer helpers
# ============================================================================================================================ #

chroms = {}
whole_genome_extracted = False

def extract_chromosome(chromosome):
    chromosome = chromosome.lower()
    if chromosome in chroms:
        system_print('Loading', chromosome, 'from cache.')
        return chroms[chromosome]
    elif chroms:
        system_print_error('Chromosome not found', chromosome)
        if whole_genome_extracted:
            return None
    c = config.Configuration()
    sequence = ''
    ref = open(c.reference)
    line = ref.readline().lower().strip()
    found = False
    while True:
        if line.startswith('>chr'):
            chrom = line[line.find('>') + 1:]
            if chrom == chromosome:
                print('extracting ' + chrom)
                while True:
                    line = ref.readline().lower().strip()
                    if line.startswith('>') or len(line) == 0:
                        print(line)
                        chroms[chromosome] = sequence
                        return sequence
                    sequence += line.upper()
        line = ref.readline().lower().strip()
        if len(line) == 0:
            break

def extract_chromosomes(chromosomes):
    c = config.Configuration()
    m = 0
    ref = open(c.reference)
    line = ref.readline().lower().strip()
    found = False
    sequence = ''
    while True:
        if line.startswith('>chr'):
            chrom = line[line.find('>') + 1:].strip().lower()
            if chrom in chromosomes:
                print('extracting ' + chrom)
                while True:
                    line = ref.readline().lower().strip()
                    if line.startswith('>') or len(line) == 0:
                        print(len(sequence), 'bases')
                        yield sequence, chrom
                        sequence = ''
                        found = True
                        m += 1
                        if m == len(chromosomes):
                            return
                        break
                    sequence += line.upper()
        # this is to avoid skipping the last line we read for the previous chromosome (header of next)
        if found:
            found = False
            continue
        line = ref.readline().lower().strip()
        if len(line) == 0:
            break

def extract_whole_genome():
    c = config.Configuration()
    global whole_genome_extracted
    if whole_genome_extracted:
        return chroms
    system_print('Extracting whole genome...')
    if not c.chromosomes:
        a = ['chr' + str(x) for x in range(1, 23)]
        a.append('chrx')
        a.append('chry')
    else:
        a = [d.lower() for d in c.chromosomes]
    for seq, chrom in extract_chromosomes(a):
        chroms[chrom] = seq
    whole_genome_extracted = True
    return chroms

