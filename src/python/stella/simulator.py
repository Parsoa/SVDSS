import io
import os
import re
import pwd
import sys
import copy
import json
import math
import time
import random
import argparse
import operator
import traceback
import subprocess

from stella import (
    bed,
    map_reduce,
    visualizer,
)

from stella.kmers import *
from stella.logger import *
from stella.chromosomes import *

print = pretty_print

import numpy as np

# ============================================================================================================================ #
# ============================================================================================================================ #
# Required arguments:
# --seed: random seed to use
# ============================================================================================================================ #
# ============================================================================================================================ #

class Simulation(map_reduce.Job):

    _name = 'TrioSimulator'
    _category = 'simulation'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = Simulation(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # Simulation
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        random.seed(c.seed + hash(self.name))
        self.tracks = self.load_structural_variations()
        self.assign_event_zygosities()
        self.export_bed(self.tracks, 'all')
        self.export_bed(self.absent, 'absent')
        self.export_bed(self.present, 'present')
        self.export_bed(self.homozygous, 'homozygous')
        self.export_bed(self.heterozygous, 'heterozygous')
        self.extract_chromosomes()
        self.index_tracks()
        self.round_robin({'chr21': self.chroms['chr21']})

    def extract_chromosomes(self):
        self.chroms = extract_whole_genome()
        self.track_index = {}
        for chrom, seq in sorted(self.chroms.items(), key = lambda x: x[0]):
            print(chrom, 'with', len(seq), 'bases')
            self.track_index[chrom] = []
            with open(os.path.join(self.get_current_job_directory(), chrom + '.fa'), 'w') as chrom_file:
                chrom_file.write('>' + chrom + '\n')
                chrom_file.write(seq.strip())
                chrom_file.write('\n')

    def load_structural_variations(self):
        c = config.Configuration()
        return self.filter_overlapping_tracks(c.tracks)

    def filter_overlapping_tracks(self, tracks):
        i = 0
        remove = []
        tracks = bed.sort_tracks(tracks)
        while i < len(tracks):
            for j in range(i + 1, len(tracks)):
                if tracks[j].chrom != tracks[i].chrom:
                    i = j
                    break
                if tracks[j].begin <= tracks[i].end:
                    remove.append(j)
                    #user_print_warning(str(tracks[j]), 'overlaps', blue(str(tracks[i])))
                    continue
                if tracks[j].begin - tracks[i].end < 1000:
                    remove.append(j)
                    #user_print_warning(str(tracks[j]), 'is too close to', blue(str(tracks[i])))
                    continue
                i = j
            if j == len(tracks) - 1:
                break
        n = 0
        for index in sorted(remove):
            tracks.pop(index - n)
            n = n + 1
        return tracks

    def assign_event_zygosities(self):
        c = config.Configuration()
        self.absent = []
        self.present = []
        self.homozygous = []
        self.heterozygous = []
        for track in self.tracks:
            m = random.randint(1, 2) - 1
            f = random.randint(1, 2) - 1
            m = 0
            f = 0
            track['allele_father'] = m
            track['allele_mother'] = f
            if m == 0 and f == 0:
                self.absent.append(track)
            else:
                self.present.append(track)
                if m == 1 and f == 1:
                    self.homozygous.append(track)
                else:
                    self.heterozygous.append(track)
        print(len(self.tracks), 'non-overlapping tracks')
        print(len(self.absent), '0|0')
        print(len(self.present), '1|1 or 1|0')
        print(len(self.homozygous), '1|1')
        print(len(self.heterozygous), '1|0')

    def index_tracks(self):
        for track in self.present:
            self.track_index[track.chrom].append(track)
        for chrom in self.track_index:
            self.track_index[chrom] = sorted(self.track_index[chrom], key = lambda x: (x.chrom, x.begin, x.end))

    def export_bed(self, tracks, name):
        print('Exporting BED file:', green(name + '.bed'))
        with open(os.path.join(self.get_current_job_directory(), name + '.bed'), 'w') as bed_file:
            for index, track in enumerate(tracks):
                if index == 0:
                    bed_file.write(track.header())
                bed_file.write(track.serialize()) 

    def transform(self, seq, chrom):
        print('simulating', chrom)
        chrom = chrom.lower()
        c = config.Configuration()
        if chrom != 'chrx' and chrom != 'chry':
            strand_1 = self.apply_events_to_chromosome(chrom, 'father')
            #strand_2 = self.apply_events_to_chromosome(chrom, 'mother')
        if chrom == 'chrx':
            print(red('Simulating chrX'))
            if self.name == 'father':
                strand_1 = self.apply_events_to_chromosome('chry', 'father')
                strand_2 = self.apply_events_to_chromosome('chrx', 'mother')
            if self.name == 'mother':
                strand_1 = self.apply_events_to_chromosome('chrx', 'father')
                strand_2 = self.apply_events_to_chromosome('chrx', 'mother')
            if self.name == 'child':
                strand_1 = self.apply_events_to_chromosome('chry', 'father')
                strand_2 = self.apply_events_to_chromosome('chrx', 'mother')
        if chrom == 'chry':
            exit()
        #self.export_diploid_chromosome_fasta(chrom, [strand_1, strand_2])
        #self.export_fastq(chrom + '_diploid')
        exit()

    def apply_events_to_chromosome(self, chrom, haplotype):
        seq = ''
        previous = 0
        for track in self.track_index[chrom.lower()]:
            seq += self.chroms[chrom][previous:track.begin]
            if track['allele_' + haplotype] == 1:
                print('applying', track, 'on', chrom) 
                if track.svtype == 'INS':
                    print('Inserting: ', track.seq.upper())
                    seq += track.seq.upper()
                if track.svtype == 'DEL':
                    pass
            else: # not present on this allele
                seq += self.chroms[chrom][track.begin: track.end]
            previous = track.end
        seq += self.chroms[chrom][previous:]
        return seq

    def export_diploid_chromosome_fasta(self, chrom, strands):
        print('Exporting FASTA file:', green(chrom + '_diploid.fa')) 
        c = config.Configuration()
        with open(os.path.join(self.get_current_job_directory(), chrom + '_diploid.fa'), 'w') as fasta_file:
            #for i in range(0, 2):
            for i in range(0, 1):
                seq = strands[i]
                #fasta_file.write('>' + chrom + '_' + str(i + 1) + '\n')
                fasta_file.write('>' + chrom + '\n')
                l = len(seq)
                k = 0
                while l >= 100:
                    fasta_file.write(seq[k: k + 100])
                    fasta_file.write('\n')
                    k += 100
                    l -= 100

    def export_fastq(self, name):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        fasta = os.path.join(self.get_current_job_directory(), name + '.fa')
        command = "pbsim --data-type CCS --depth 30 --sample-fastq /share/hormozdiarilab/Codes/Fulminata/data/Samples/HG00731/HG00731_20190925_EEE_m54329U_190528_231241.Q20.fastq --prefix {} --length-max 20000 --length-min 2000 {}".format(name, fasta)
        print(command)
        if self.name == 'child':
            output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT, cwd = self.get_current_job_directory())

    def reduce(self):
        c = config.Configuration()
        return self.present

    def get_current_job_directory(self):
        return os.path.abspath(os.path.join(self.get_output_directory(), self._name, self.name))

# ============================================================================================================================ #
# ============================================================================================================================ #
# Required arguments:
# --seed: random seed to use
# ============================================================================================================================ #
# ============================================================================================================================ #

class ChildSimulation(Simulation):

    _name = 'TrioSimulator'
    _category = 'simulation'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = ChildSimulation(**kwargs)
        job.execute()

    # ============================================================================================================================ #
    # Simulation
    # ============================================================================================================================ #

    def load_inputs(self):
        c = config.Configuration()
        self.export_bed(self.present, 'present')
        self.extract_chromosomes()
        self.index_tracks()
        self.round_robin({'chr21': self.chroms['chr21']})

    def apply_events_to_chromosome(self, chrom, haplotype):
        seq = self.chroms[chrom][self.track.begin - 10000: self.track.begin]
        if self.track.svtype == 'INS':
            seq += self.track.seq
        if self.track.svtype == 'DEL':
            track.end = track.begin + abs(int(track.svlen))
            pass
        seq += self.chroms[chrom][self.track.end: self.track.end + 10000]
        return seq

# ============================================================================================================================ #
# ============================================================================================================================ #
# Required arguments:
# --seed: random seed to use
# ============================================================================================================================ #
# ============================================================================================================================ #

class LociSimulator(Simulation):

    def load_inputs(self):
        c = config.Configuration()
        self.export_bed(self.present, 'present')
        self.extract_chromosomes()
        self.round_robin({'chr21': self.chroms['chr21']})

    def apply_events_to_chromosome(self, chrom, haplotype):
        seq = ''
        previous = 0
        seq += self.chroms[chrom][previous:track.begin]
        previous = track.end
        if track['allele_' + haplotype] == 1:
            print('applying ', track, ' to chr21')
            track.begin = len(seq)
            if track.svtype == 'INS':
                seq += track.seq
                track.end = track.begin + 1
            if track.svtype == 'DEL':
                track.end = track.begin + abs(int(track.svlen))
                pass
        else: # not present on this allele
            track.begin = len(seq)
            print('not applying ', track, ' to chr21')
            seq += self.chroms[chrom][track.begin: track.end]
            if track.svtype == 'INS':
                track.end = track.begin + 1
                pass
            if track.svtype == 'DEL':
                track.end = track.begin + abs(int(track.svlen))
                pass
        seq += self.chroms[chrom][previous:]
        self.export_bed(self.present, 'actual')
        return seq

# ============================================================================================================================ #
# ============================================================================================================================ #
# Required arguments:
# --seed: random seed to use
# ============================================================================================================================ #
# ============================================================================================================================ #

class TrioSimulator(map_reduce.Job):

    _name = 'TrioSimulator'
    _category = 'simulation'
    _previous_job = None

    # ============================================================================================================================ #
    # Launcher
    # ============================================================================================================================ #

    @staticmethod
    def launch(**kwargs):
        job = TrioSimulator(**kwargs)
        job.execute()

    def load_inputs(self):
        c = config.Configuration()
        extract_whole_genome()
        self.tracks = {'father': True, 'mother': True}
        self.round_robin(self.tracks)
        self.reduce()

    def transform(self, track, track_name):
        job = Simulation(name = track_name)
        tracks = job.execute()
        path = os.path.join(self.get_current_job_directory(), track_name)
        with open(os.path.join(path, 'present.bed'), 'w') as bed_file:
            for i, track in enumerate(tracks):
                if i == 0:
                    bed_file.write(track.header())
                bed_file.write(track.serialize())
        return None

    def reduce(self):
        tracks = []
        _tracks = {}
        for parent in ['father', 'mother']:
            path = os.path.join(self.get_current_job_directory(), parent, 'all.bed')
            _tracks[parent] = bed.load_tracks_from_file(path)
        random.seed(92106429)
        for f, m in zip(_tracks['father'], _tracks['mother']):
            if f.chrom != 'chr21':
                continue
            t = copy.deepcopy(f)
            # chrX is always inherited from mother
            if f.chrom.lower() == 'chrx':
                t['allele_father'] = '.'
                t['allele_mother'] = random.choice([int(m['allele_father']), int(m['allele_mother'])])
            # chrY always inherited from father and has a single haplotype
            elif f.chrom.lower() == 'chry':
                t['allele_father'] = int(f['allele_father'])
                t['allele_mother'] = '.'
            # other chromosomes are normal
            else:
                r = random.choice([0, 1])
                t['allele_father'] = r
                t['allele_mother'] = r
                #t['allele_father'] = random.choice([int(f['allele_father']), int(f['allele_mother'])])
                #t['allele_mother'] = random.choice([int(m['allele_father']), int(m['allele_mother'])])
                if r == 1:
                    tracks.append(t)
        job = ChildSimulation(name = 'child', present = tracks)
        tracks = job.execute()
        exit()

