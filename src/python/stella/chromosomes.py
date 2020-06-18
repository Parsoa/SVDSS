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
import argparse
import traceback

from stella import (
    config,
)

from stella.logger import *

# ============================================================================================================================ #
# kmer helpers
# ============================================================================================================================ #

chroms = {}
whole_genome_extracted = False

def extract_chromosomes(chromosomes):
    c = config.Configuration()
    m = 0
    ref = open(c.reference)
    line = ref.readline().strip()
    found = False
    sequence = ''
    while True:
        if line.startswith('>chr'):
            chrom = line[line.find('>') + 1:].strip()
            if chrom in chromosomes:
                while True:
                    line = ref.readline().strip()
                    if line.startswith('>') or len(line) == 0:
                        system_print('Extracted ' + chrom + ' with', len(sequence), 'bases.')
                        yield sequence, chrom
                        sequence = ''
                        found = True
                        m += 1
                        if m == len(chromosomes):
                            return
                        break
                    sequence += line.upper()
        if found:
            found = False
            continue
        line = ref.readline().strip()
        if len(line) == 0:
            break

def extract_whole_genome():
    c = config.Configuration()
    if c.reduce:
        a = ['chr' + str(x) for x in range(1, 23)]
        a.append('chrX')
        a.append('chrY')
        for chrom in a:
            chroms[chrom] = 'ATCG'
        return chroms
    global whole_genome_extracted
    if whole_genome_extracted:
        return chroms
    system_print('Loading reference assembly..')
    if not c.chromosomes:
        a = ['chr' + str(x) for x in range(1, 23)]
        a.append('chrX')
        a.append('chrY')
    else:
        a = [d for d in c.chromosomes]
    print(a)
    for seq, chrom in extract_chromosomes(a):
        chroms[chrom] = seq
    whole_genome_extracted = True
    return chroms

