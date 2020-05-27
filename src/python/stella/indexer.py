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

from suffix_tree import Tree

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class SuffixTreeIndexCreator(map_reduce.SingleCoreJob):

    _name = 'SuffixTreeIndexCreator'
    _category = 'preprocessing'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = SuffixTreeIndexCreator(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        cpp_dir = os.path.join(os.path.dirname(__file__), '../../cpp')
        command = os.path.join(cpp_dir, "counter.out") + " " + c.father + " " + c.mother + " " + c.child + " " + self.get_current_job_directory() + " " + "chr20"
        print(command)
        output = subprocess.call(command, shell = True)
        exit()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class SuffixTreePostProcessor(map_reduce.SingleCoreJob):

    _name = 'SuffixTreePostProcessor'
    _category = 'preprocessing'
    _previous_job = SuffixTreeIndexCreator

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = SuffixTreePostProcessor(**kwargs)
        job.execute()

    def decode_base(self, base):
        if base == 65:
            return 'A'
        if base == 67:
            return 'C'
        if base == 71:
            return 'G'
        if base == 84:
            return 'T'
        return 'N'

    def load_inputs(self):
        reads = []
        i = 0
        while True:
            batch = self.load_previous_output_batch(i)
            if not batch:
                break
            print(cyan('Reading batch', i))
            i += 1
            for read in batch:
                s = "".join([self.decode_base(int(x)) for x in read.split()])
                reads.append(s)
        reads = sorted(reads, key = lambda x: len(x), reverse = True)
        print('Merged and sorted all reads:', len(reads))
        self.export_reads_as_fasta(reads)
        with open(os.path.join(self.get_current_job_directory(), 'reads.json'), 'w') as json_file:
            json.dump(reads, json_file, indent = 4)
        _reads = {i: r for i, r in enumerate(reads)}
        print('Building suffix tree..')
        tree = Tree(_reads)
        print('Done..')
        remove = {}
        n = 0
        while True:
            print('Iteration', n, '.', len(_reads), 'sequences remaining.')
            l = len(_reads)
            for i in range(len(reads)):
                if not i in remove:
                    for _id, path in tree.find_all(reads[i]):
                        if _id != i:
                            remove[i] = True
                            break
            for i in remove:
                if i in _reads:
                    _reads.pop(i, None)
            if len(_reads) == l:
                break
            l = len(_reads)
        print(len(reads), len(remove))

    def export_reads_as_fasta(self, reads):
        with open(os.path.join(self.get_current_job_directory(), 'reads_unfiltered.fa'), 'w') as fasta_file:
            for i, read in enumerate(reads):
                fasta_file.write('>read' + str(i) + '\n')
                fasta_file.write(read + '\n')
        print('Dumped reads.')

    def reduce(self):
        with open(os.path.join(self.get_current_job_directory(), 'reads.json'), 'r') as json_file:
            reads = json.load(json_file)
            visualizer.histogram([[len(x) for x in reads]], 'Novel Sequence Length Distribution', self.get_current_job_directory(), 'Lenght of Sequence', 'Number if Sequences')

