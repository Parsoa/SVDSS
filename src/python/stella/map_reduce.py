from __future__ import print_function

import io
import gc
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import atexit
import argparse
import traceback

from shutil import copyfile

from stella import (
    bed,
    config,
    counttable,
)

from stella.debug import *
from stella.kmers import *
from stella.logger import *
print = pretty_print

def on_exit(job):
    print(green('job', job.index, 'exiting'))

# ============================================================================================================================ #
# Job class, describes a MapReducer job
# Will apply a transformation to a library of structural variation:
# 1. Divides the library into a number of `batches`
# 2. Each batch includes one or more `tracks` (each track is a structural variation)
# 3. Applies the function `transform` to each track and outputs the result
# 4. Merges the transformed tracks into a whole
# ============================================================================================================================ #

class Job(object):

    def __init__(self, **kwargs):
        c = config.Configuration()
        self.index = -1
        self.batch = {}
        self.children = {}
        self.resume_from_reduce = c.reduce
        for k, v in kwargs.items():
            setattr(self, k, v)

    def execute(self):
        start = time.clock() 
        c = config.Configuration()
        self.create_output_directories()
        self.load_inputs()
        if not self.resume_from_reduce: 
            self.distribute_workload()
            self.wait_for_children()
        else:
            system_print_high('Resuming from reduce...')
        output = self.reduce()
        end = time.clock()
        system_print_high('Stage ' + self._name + ' finished. Execution time', end - start)
        return output

    def load_inputs(self):
        tracks = self.load_previous_job_results()
        self.round_robin(tracks)

    def load_previous_job_results(self, job = None):
        if not job:
            path = os.path.join(self.get_previous_job_directory(), 'batch_merge.json')
        else:
            path = os.path.join(os.path.join(self.get_current_job_directory(), '..', job._name), 'batch_merge.json')
        with open(path, 'r') as json_file:
            return json.load(json_file)

    def round_robin(self, tracks, filter_func = lambda x: False):
        c = config.Configuration()
        n = 0
        self.batch = {}
        self.num_threads = 0
        for track in tracks:
            if filter_func(tracks[track]):
                continue
            index = n % c.threads
            if not index in self.batch:
                self.batch[index] = {}
            self.batch[index][track] = tracks[track]
            system_print('assigned ', track, ' to ', index)
            n = n + 1
            self.num_threads = min(c.threads, n)

    def distribute_workload(self):
        for index in range(0, self.num_threads):
            pid = os.fork()
            if pid == 0:
                self.index = index
                self.run_batch(self.batch[index])
                exit()
            else:
                self.children[pid] = index
                system_print('spawned child', '{:2d}'.format(index), ':', pid)
        system_print('done distributing workload')

    def run_batch(self, batch):
        c = config.Configuration()
        remove = {}
        n = 0
        start = time.time()
        for track in batch:
            try:
                batch[track] = self.transform(batch[track], track)
                if batch[track] == None:
                    remove[track] = True
            except Exception as e:
                system_print_error(track)
                system_print_error(traceback.format_exc())
                remove[track] = True
                debug_breakpoint()
            n = n + 1
            t = time.time()
            p = float(n) / len(batch)
            eta = (1.0 - p) * ((1.0 / p) * (t - start)) / 3600
            system_print('{:2d}'.format(self.index), 'progress:', '{:7.5f}'.format(p), 'ETA:', '{:8.6f}'.format(eta))
            if n % 1000 == 0:
                gc.collect()
        for track in remove:
            batch.pop(track, None)
        # if there is no output, don't write anything
        self.output_batch(batch)
        self.on_exit_worker()
        exit()

    def on_exit_worker(self):
        pass

    def transform(self, track, track_name):
        return track

    def output_batch(self, batch):
        json_file = open(os.path.join(self.get_current_job_directory(), 'batch_' + str(self.index) + '.json'), 'w')
        json.dump(batch, json_file, sort_keys = True, indent = 4)
        json_file.close()

    def wait_for_children(self):
        while self.children:
            (pid, e) = os.wait()
            index = self.children[pid]
            self.children.pop(pid, None)
            if os.path.isfile(os.path.join(self.get_current_job_directory(), 'batch_' + str(index) + '.json')):
                system_print_high('pid', '{:5d}'.format(pid) + ', index', '{:2d}'.format(index), 'finished,', '{:2d}'.format(len(self.children)), 'remaining')
            else:
                system_print_high('pid', '{:5d}'.format(pid) + ', index', '{:2d}'.format(index), 'finished didn\'t produce output,', len(self.children), 'remaining')
        system_print_high('All forks done, merging output...')

    def reduce(self):
        c = config.Configuration()
        output = {}
        for i in range(0, self.num_threads):
            path = os.path.join(self.get_current_job_directory(), 'batch_' + str(i) + '.json')
            if os.path.isfile(path):
                with open(path, 'r') as json_file:
                    batch = json.load(json_file)
                    output.update(batch)
        with open(os.path.join(self.get_current_job_directory(), 'batch_merge.json'), 'w') as json_file:
            json.dump(output, json_file, sort_keys = True, indent = 4)
        return output

    def load_output(self):
        for i in range(0, self.num_threads):
            system_print('reading batch', i)
            yield self.load_output_batch(i)

    def load_output_batch(self, index):
        path = os.path.join(self.get_current_job_directory(), 'batch_' + str(index) + '.json')
        if not os.path.isfile(path):
            system_print_error('didn\'t find batch', index)
            return {}
        with open(path, 'r') as json_file:
            output = json.load(json_file)
            return output

    # ============================================================================================================================ #
    # misc helpers
    # ============================================================================================================================ #

    def load_tracks(self, name = 'all.bed'):
        c = config.Configuration()
        if c.simulation:
            return bed.load_tracks_from_file_as_dict(os.path.join(self.get_simulation_directory(), name))
        else:
            return bed.load_tracks_from_file_as_dict(c.bed)

    def load_reference_counts_provider(self):
        c = config.Configuration()
        self.reference_counts_provider = counttable.JellyfishCountsProvider(c.jellyfish)

    def unload_reference_counts_provider(self):
        del self.reference_counts_provider

    # ============================================================================================================================ #
    # filesystem helpers
    # ============================================================================================================================ #

    def get_output_directory(self):
        c = config.Configuration()
        return os.path.abspath(c.workdir)

    def get_current_job_directory(self):
        return os.path.abspath(os.path.join(self.get_output_directory(), self._name))

    def get_previous_job_directory(self):
        c = config.Configuration()
        if self._previous_job:
            return os.path.abspath(os.path.join(self.get_output_directory(), self._previous_job._name))
        else:
            return None

    def create_output_directories(self):
        path = self.get_output_directory()
        if not os.path.exists(path):
            os.makedirs(path)
        path = self.get_current_job_directory()
        if not os.path.exists(path):
            os.makedirs(path)

# ============================================================================================================================ #
# ============================================================================================================================ #
# Base class for every job that is a direct part of the genotyping process
# ============================================================================================================================ #
# ============================================================================================================================ #

class GenomeDependentJob(Job):

    def get_current_job_directory(self):
        c = config.Configuration()
        if c.simulation:
            s = Job.get_current_job_directory(self)
            print(yellow(s))
            return s
        else:
            return os.path.abspath(os.path.join(self.get_output_directory(), self._name))

# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #
# ============================================================================================================================ #

class TrackExportHelper(Job):

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    _name = 'TrackExportHelper'
    _category = 'genotyping'
    _previous_job = None

    @staticmethod
    def launch(**kwargs):
        job = TrackExportHelper(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # MapReduce overrides
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.resume_from_reduce = False
        self.round_robin(self.tracks)

    def transform(self, track, track_name):
        with open(os.path.join(self.get_current_job_directory(), track_name + '.json'), 'w') as json_file:
            json.dump(track, json_file, indent = 4)
        return None

    def get_current_job_directory(self):
        return self.current_job_directory

    def output_batch(self, batch):
        pass

    def reduce(self):
        pass
