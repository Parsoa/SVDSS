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

from MUSCython import MultiStringBWTCython as MSBWT

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class BwtDiffJob(map_reduce.Job):

    _name = 'BwtDiffJob'
    _category = 'preprocessing'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = BwtDiffJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def round_robin(self):
        c = config.Configuration()
        self.num_threads = 1 if c.debug else 8
        for index in range(self.num_threads):
            self.batch[index] = {}

    def run_batch(self, batch):
        self.transform()

    def load_inputs(self):
        c = config.Configuration()
        print('Loading BWT ' + c.bwt)
        self.msbwt = MSBWT.loadBWT(c.bwt)
        self.round_robin()

    def transform(self):
        c = config.Configuration()
        self.mismatched_strings = {}
        #self.bam_file = pysam.AlignmentFile(c.bam, "rb")
        #print('Thread', self.index, 'beginning..')
        #n = 0
        #reads = self.bam_file.fetch(until_eof = True)# contig = 'chr21')
        #for read in reads:
        #    if n % self.num_threads == self.index:
        #        if read.query_sequence:
        #            self.ping_pong_search(read.query_sequence)
        #    n += 1
        #    if n % 100 == 0:
        #        print('\rThread {:>4}, Progress {:>10d}. Found {:>10d}.'.format(self.index, n, len(self.mismatched_strings)))
        with open(c.fastq[0]) as fastq_file:
            n = 0
            m = 0
            line = fastq_file.readline()
            while line:
                if n % 4 == 1:
                    if m % self.num_threads == self.index:
                        self.ping_pong_v2(line)
                    m += 1
                line = fastq_file.readline()
                n += 1
                if n % 1000 == 0:
                    print("\rThread", self.index, ". Progress:", n, '. Found', len(self.mismatched_strings))
        self.output_batch(self.mismatched_strings)
    
    def ping_pong_v2(self, line):
        line = line.strip()
        l = len(line)
        begin = l - 1
        while (begin >= 0):
            # Backward search
            tLow, tHigh = self.msbwt.findIndicesOfStr(line[begin])
            debug_log("BS from " + line[begin] + " (" + str(begin) + "): [" + str(tHigh) + ", " + str(tLow) +  "," + str(tHigh - tLow) + "]")
            while tHigh - tLow != 0 and begin > 0:
                begin -= 1
                tLow, tHigh = self.msbwt.findIndicesOfStr(line[begin], (tLow, tHigh))
                debug_log("- BE with " + line[begin] + " (" + str(begin) + "): [" + str(tHigh) + ", " + str(tLow) + ", " + str(tHigh - tLow) + "]")
            if begin == 0 and tHigh - tLow != 0:
                break
            debug_log("Mismatch " + line[begin] + " (" + str(begin) + ")")
            # Forward search
            end = begin
            debug_log("FS from " + line[end] + " (" + str(end) + ")")
            m = self.msbwt.countOccurrencesOfSeq(line[begin: end + 1])
            while m != 0:
                end += 1
                m = self.msbwt.countOccurrencesOfSeq(line[begin: end + 1])
                debug_log("- FE from " + line[end] + " (" + str(end) + ")")
            s = line[begin: end + 1]
            assert self.msbwt.countOccurrencesOfSeq(s) == 0
            if not s in self.mismatched_strings:
                self.mismatched_strings[s] = 0
            self.mismatched_strings[s] += 1
            debug_sleep(1)

    def ping_pong_v1(self, line):
        line = line.strip()
        l = len(line)
        i = l - 1
        tLow, tHigh = self.msbwt.findIndicesOfStr(line[i])
        i -= 1
        while i >= 0:
            tLow, tHigh = self.msbwt.findIndicesOfStr(line[i], (tLow, tHigh))
            if tHigh - tLow == 0:
                break
            i -= 1
        if i != -1:
            end = l - 1
            begin = i 
            prev_end = end
            while begin >= 0:
                debug_log("== begin: " + str(begin) + ", end: " + str(end))
                # Forward search
                # find smallest prefix of [begin, end] that is still a mismatch
                while end > begin:
                    debug_log("-E interval [" + str(begin) + ", " + str(end) + "]")
                    m = self.msbwt.countOccurrencesOfSeq(line[begin: end + 1])
                    if m == 0:
                        prev_end = end
                    else:
                        debug_log("## fixed end at " + str(prev_end))
                        s = line[begin: prev_end + 1]
                        assert self.msbwt.countOccurrencesOfSeq(s) == 0
                        if s not in self.mismatched_strings:
                            self.mismatched_strings[s] = 0
                        self.mismatched_strings[s] += 1
                        break
                    end -= 1
                begin -= 1
                # find the next mismatch and repeat
                s = False
                while begin >= 0:
                    debug_log("-B interval [" + str(begin) + ", " + str(end) + "]")
                    tLow, tHigh = self.msbwt.findIndicesOfStr(line[begin], (tLow, tHigh)) if s else self.msbwt.findIndicesOfStr(line[begin: end + 1])
                    debug_log("%%", tHigh - tLow)
                    s = True
                    if tHigh - tLow == 0:
                        assert self.msbwt.countOccurrencesOfSeq(line[begin: end + 1]) == 0
                        debug_log("@@ fixed begin at " + str(begin))
                        break
                    begin -= 1
                debug_sleep(1)

    def reduce(self):
        pass

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class BwtPostProcessor(map_reduce.SingleCoreJob):

    _name = 'BwtPostProcessor'
    _category = 'preprocessing'
    _previous_job = BwtDiffJob

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = BwtPostProcessor(**kwargs)
        job.execute()

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

