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
import argparse
import operator
import traceback
import subprocess

from threading import Thread, Lock

from stella import (
    bed,
    config,
    map_reduce,
    visualizer,
)

from stella.debug import *
from stella.kmers import *
from stella.logger import *
from stella.chromosomes import *

print = pretty_print

import pysam

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class FmdIndexCreator(map_reduce.Job):

    _name = 'FmdIndexCreator'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = FmdIndexCreator(**kwargs)
        job.execute()

    def load_inputs(self):
        c = config.Configuration()
        exec_path = '../../cpp/fmd/main.out'
        command = ' '.join([exec_path, 'index', c.fastq, '>', c.fmd])

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class FmdShortSequenceScanner(map_reduce.Job):

    _name = 'FmdShortSequenceScanner'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = FmdIndexCreator(**kwargs)
        job.execute()

    def load_inputs(self):
        c = config.Configuration()
        chroms = extract_whole_genome()
        self.load_sequences()
        print('Scanning for', len(self.kmers), 'kmers. Min length:', self.min_length, 'Max length:', self.max_length)
        self.transform(chroms[c.chromosomes[0]], c.chromosomes[0])
        self.reduce()
        exit()

    def read_fastq_entry(self, fastq_file):
        fastq = []
        line = fastq_file.readline()
        if not line:
            return None
        fastq.append(line)
        fastq.append(fastq_file.readline())
        fastq.append(fastq_file.readline())
        fastq.append(fastq_file.readline())
        return fastq

    def load_sequences(self):
        c = config.Configuration()
        self.kmers = {}
        self.max_length = 0
        self.min_length = 1000
        with open(c.fastq[0]) as fastq_file:
            n = 0
            fastq = self.read_fastq_entry(fastq_file)
            while fastq:
                kmer = canonicalize(fastq[1].strip())
                self.kmers[kmer] = {'header': fastq[0], 'loci': [], 'seq': fastq[1].strip()}
                self.min_length = len(kmer) if self.min_length > len(kmer) else self.min_length
                self.max_length = len(kmer) if self.max_length < len(kmer) else self.max_length
                n += 1
                fastq = self.read_fastq_entry(fastq_file)
                if n % 1000000 == 0:
                    print('Processed', n, 'kmers..')

    def transform(self, sequence, chrom):
        c = config.Configuration()
        t = time.time()
        l = len(sequence)
        for index in range(0, l - 100):
            for k in range(10, 40):
                kmer = canonicalize(sequence[index: index + k])
                if kmer in self.kmers:
                    self.kmers[kmer]['loci'].append(index)
            if index % 1000000 == 1:
                s = time.time()
                p = index / float(len(sequence))
                e = (1.0 - p) * (((1.0 / p) * (s - t)) / 3600)
                system_print('{:5}'.format(chrom), 'progress:', '{:12.10f}'.format(p), 'took:', '{:14.10f}'.format(s - t), 'ETA:', '{:12.10f}'.format(e))

    def reduce(self):
        n = 0
        with open(os.path.join(self.get_current_job_directory(), 'solution.short.fastq'), 'w') as fastq_file:
            for kmer in self.kmers:
                if self.kmers[kmer]['loci']:
                    n += 1
                for pos in self.kmers[kmer]['loci']:
                    fastq_file.write('\t'.join(['chr21', str(pos), self.kmers[kmer]['seq'], self.kmers[kmer]['header']]) + '\n')
        print('Mapped', n, 'sequences.')
