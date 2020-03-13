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
    config,
    map_reduce,
)

from stella.debug import *
from stella.kmers import *
from stella.logger import *
from stella.chromosomes import *

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class BamProcessingJob(map_reduce.Job):

    _name = 'ReferenceKmerExtractor'
    _category = 'preprocessing'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.tracks = c.tracks
        if not c.reduce:
            extract_whole_genome()
        self.round_robin(self.tracks)
        self.load_reference_counts_provider()
