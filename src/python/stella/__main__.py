from stella import (
    config,
    indexer,
    simulator,
    map_reduce,
    preprocessor,
)

from stella.bwt import *
from stella.debug import *
from stella.kmers import *
from stella.logger import *
from stella.chromosomes import *

# ============================================================================================================================ #

def load_tracks(filter_overlap = True):
    job = preprocessor.TrackPreprocessorJob(filter_overlap = filter_overlap)
    tracks = job.execute()
    config.Configuration.update({'tracks': tracks})

def fmd_index():
    job = FmdIndexCreator()
    job.execute()

def bwt():
    job = BwtDiffJob()
    job.execute()
    job = BwtPostProcessor()
    #job.execute()

def index():
    job = indexer.SuffixTreeIndexCreator()
    job.execute()
    job = indexer.SuffixTreePostProcessor()
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
    if c.command == 'index':
        if c.bwt:
            bwt()
        if c.fmd:
            fmd()
        else:
            index()
    if c.command == 'simulate':
        simulate()
