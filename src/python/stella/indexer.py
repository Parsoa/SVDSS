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
)

from stella.debug import *
from stella.kmers import *
from stella.logger import *
from stella.chromosomes import *
print = pretty_print

import pysam

import numpy as np

from suffix_trees import STree

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class SuffixTreeIndexCreator(map_reduce.Job):

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
        self.alignments = pysam.AlignmentFile(c.bam, "rb", check_sq=False)
        seqs = []
        n = 0
        for read in self.alignments.fetch(until_eof = True):
            seqs.append(read.query_sequence)
            n += 1
            if n == 1000:
                break
        print('Creating tree on', len(seqs), 'sequences..')
        tree = STree.STree(seqs)
        print('Done')
        debug_breakpoint()
        exit()

