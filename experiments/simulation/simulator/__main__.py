from stella import (
    config,
    simulator,
    map_reduce,
    preprocessor,
)

from stella.fmd import *
from stella.debug import *
from stella.kmers import *
from stella.logger import *
from stella.chromosomes import *

# ============================================================================================================================ #

def load_tracks(filter_overlap = True):
    job = preprocessor.TrackPreprocessorJob(filter_overlap = filter_overlap)
    tracks = job.execute()
    config.Configuration.update({'tracks': tracks})

def scan():
    job = FmdShortSequenceScanner()
    job.execute()

def fmd_index():
    job = FmdIndexCreator()
    job.execute()

def simulate():
    load_tracks()
    job = simulator.TrioSimulator()
    job.execute()

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# Main
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

if __name__ == '__main__':
    config.init()
    c = config.Configuration()
    if c.command == 'scan':
        scan()
    if c.command == 'index':
        fmd_index()
    if c.command == 'search':
        fmd_index()
    if c.command == 'simulate':
        simulate()
