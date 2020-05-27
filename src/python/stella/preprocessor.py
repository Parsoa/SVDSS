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

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class TrackPreprocessorJob(map_reduce.Job):

    _name = 'TrackPreprocessorJob'
    _category = 'preprocessing'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = TrackPreprocessorJob(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    SUPPORTED_SVTYPES = ['DEL', 'INS', 'INV', 'ALU', 'MEI']

    def execute(self):
        c = config.Configuration()
        tracks = {}
        for path in c.bed:
            tracks.update(bed.load_tracks_from_file_as_dict(path, parse_header = True))
        print('Loaded', len(tracks), 'tracks.')
        tracks = bed.sort_tracks(tracks)
        svtypes = {track.svtype: True for track in tracks}
        merged_tracks = []
        for svtype in svtypes:
            if svtype in TrackPreprocessorJob.SUPPORTED_SVTYPES:
                _tracks = [track for track in tracks if track.svtype == svtype]
                print('Merging', svtype + '.', len(_tracks), 'tracks.')
                _tracks = bed.filter_overlapping_tracks(_tracks, svtype)
                print('Merged.', len(_tracks), 'tracks remaining.')
                merged_tracks += _tracks
        tracks = bed.sort_tracks(merged_tracks)
        tracks = {track.id: track for track in tracks if track.svtype in TrackPreprocessorJob.SUPPORTED_SVTYPES}
        print('Removed overlapping tracks.', len(tracks), 'non-overlapping tracks.')
        return tracks

