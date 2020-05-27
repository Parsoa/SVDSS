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

