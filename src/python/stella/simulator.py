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
        self.round_robin({chrom: self.chroms[chrom] for chrom in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']})

    def extract_chromosomes(self):
        self.chroms = extract_whole_genome()

    def load_structural_variations(self):
        c = config.Configuration()
        return bed.sort_tracks(c.tracks)

    def assign_event_zygosities(self):
        c = config.Configuration()
        self.absent = []
        self.present = []
        self.homozygous = []
        self.heterozygous = []
        for track in self.tracks:
            f = random.randint(1, 2) - 1
            m = random.randint(1, 2) - 1
            track['allele_father'] = str(f)
            track['allele_mother'] = str(m)
            if m == 0 and f == 0:
                self.absent.append(track)
                track['genotype'] = '0/0'
            else:
                self.present.append(track)
                if m == 1 and f == 1:
                    track['genotype'] = '1/1'
                    self.homozygous.append(track)
                else:
                    track['genotype'] = '1/0'
                    self.heterozygous.append(track)
        print(len(self.tracks), 'non-overlapping tracks')
        print(len(self.absent), '0|0')
        print(len(self.present), '1|1 or 1|0')
        print(len(self.homozygous), '1|1')
        print(len(self.heterozygous), '1|0')

    def index_tracks(self):
        self.track_index = {}
        for track in self.present:
            if not track.chrom in self.track_index:
                self.track_index[track.chrom] = []
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
            strand_2 = self.apply_events_to_chromosome(chrom, 'mother')
        if chrom == 'chrx':
            print(red('Simulating chrX'))
            if self.name == 'father':
                strand_1 = self.apply_events_to_chromosome('chry', 'father')
                strand_2 = self.apply_events_to_chromosome('chrx', 'mother')
            if self.name == 'mother':
                strand_1 = self.apply_events_to_chromosome('chrx', 'father')
                strand_2 = self.apply_events_to_chromosome('chrx', 'mother')
            # child is female
            if self.name == 'child':
                strand_1 = self.apply_events_to_chromosome('chrx', 'father')
                strand_2 = self.apply_events_to_chromosome('chrx', 'mother')
        if chrom == 'chry':
            exit()
        self.export_diploid_chromosome_fasta(chrom, [strand_1, strand_2])
        self.export_fastq(chrom + '_diploid_1')
        self.export_fastq(chrom + '_diploid_2')
        exit()

    def apply_events_to_chromosome(self, chrom, haplotype):
        seq = ''
        previous = 0
        for track in self.track_index[chrom.lower()]:
            seq += self.chroms[chrom][previous:track.begin]
            if track['allele_' + haplotype] == '1':
                print('applying', track, 'on', chrom) 
                if track.svtype == 'INV':
                    seq += reverse_complement(self.chroms[chrom][track.begin: track.end])
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
        for i in range(0, 2):
            with open(os.path.join(self.get_current_job_directory(), chrom + '_diploid_' + str(i + 1) + '.fa'), 'w') as fasta_file:
                fasta_file.write('>' + chrom + '_' + str(i + 1) + '\n')
                seq = strands[i]
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
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT, cwd = self.get_current_job_directory())

    def reduce(self):
        c = config.Configuration()
        return self.present

    def get_current_job_directory(self):
        return os.path.abspath(os.path.join(self.get_output_directory(), self.name))

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

    # ============================================================================================================================ ## Launcher
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
        self.assign_event_zygosities()
        self.export_bed(self.tracks, 'all')
        self.export_bed(self.absent, 'absent')
        self.export_bed(self.present, 'present')
        self.export_bed(self.homozygous, 'homozygous')
        self.export_bed(self.heterozygous, 'heterozygous')
        self.extract_chromosomes()
        self.index_tracks()
        self.round_robin({chrom: self.chroms[chrom] for chrom in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']})

    def assign_event_zygosities(self):
        c = config.Configuration()
        self.absent = []
        self.present = []
        self.homozygous = []
        self.heterozygous = []
        for track in self.tracks:
            f = int(track['allele_father'])
            m = int(track['allele_mother'])
            if m == 0 and f == 0:
                self.absent.append(track)
                track['genotype'] = '0/0'
            else:
                self.present.append(track)
                if m == 1 and f == 1:
                    track['genotype'] = '1/1'
                    self.homozygous.append(track)
                else:
                    track['genotype'] = '1/0'
                    self.heterozygous.append(track)
        print(len(self.tracks), 'non-overlapping tracks')
        print(len(self.absent), '0|0')
        print(len(self.present), '1|1 or 1|0')
        print(len(self.homozygous), '1|1')
        print(len(self.heterozygous), '1|0')

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
        #########extract_whole_genome()
        self.tracks = {'father': True, 'mother': True}
        self.round_robin(self.tracks)

    def transform(self, track, track_name):
        job = Simulation(name = track_name)
        tracks = job.execute()
        return None

    def reduce(self):
        c = config.Configuration()
        tracks = []
        _tracks = {}
        chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']
        for parent in ['father', 'mother']:
            path = os.path.join(self.get_current_job_directory(), parent, 'all.bed')
            _tracks[parent] = bed.load_tracks_from_file(path)
        random.seed(92106429)
        haploids = ['allele_mother', 'allele_father']
        mother_haploid = random.choice([0, 1])
        father_haploid = random.choice([0, 1])
        tracks_per_chrom = {}
        for track in _tracks['father']:
            if track.chrom not in tracks_per_chrom:
                tracks_per_chrom[track.chrom] = 0
            tracks_per_chrom[track.chrom] += 1
        json_print(tracks_per_chrom)
        for f, m in zip(_tracks['father'], _tracks['mother']):
            t = copy.deepcopy(f)
            # Different headers cause problem when writing BED files
            if t.svlen == 'INV':
                b = bed.BedTrack(chrom = t.chrom, begin = t.begin, end = t.end) 
                b['svtype'] = 'INV'
                b['svlen'] = int(t.svtype)
                b['seq'] = '.'
                b['svclass'] = 'Inversion'
                b['is_trf'] = 'False'
                b['callset'] = 'pacbio'
                b['genotype'] = t.seq
                b['allele_father'] = t.svclass
                b['allele_mother'] = t.is_trf
                m.allele_father = m.svclass
                m.allele_mother = m.is_trf
                f.allele_father = f.svclass
                f.allele_mother = f.is_trf
                t = b
            if t.chrom not in chromosomes:
                continue
            if t.id not in c.tracks:
                print('Skipping', t.id)
                continue
            #print(t['allele_mother'], t['allele_father'])
            t['allele_mother'] = m[haploids[mother_haploid]]
            t['allele_father'] = f[haploids[father_haploid]]
            t['mother_haploid'] = mother_haploid
            t['father_haploid'] = father_haploid
            if t['allele_father'] == '0' and t['allele_mother'] == '0':
                r = random.choice([0, 2])
                if r == 0:
                    t['allele_father'] = '0'
                    t['allele_mother'] = '1'
                if r == 1:
                    t['allele_father'] = '1'
                    t['allele_mother'] = '0'
                if r == 2:
                    t['allele_father'] = '1'
                    t['allele_mother'] = '1'
            tracks.append(t)
            # recombination
            choices = [1, 2]
            p = 10.0 / tracks_per_chrom[t.chrom]
            weights = [p, 1.0 - p]
            r = np.random.choice(choices, p = weights)
            if r == 1:
                print('Changing mother haploid')
                mother_haploid = (mother_haploid + 1) % 2
            r = np.random.choice(choices, p = weights)
            if r == 1:
                print('Changing father haploid')
                father_haploid = (father_haploid + 1) % 2
        job = ChildSimulation(name = 'child', tracks = tracks)
        tracks = job.execute()

    def get_current_job_directory(self):
        return os.path.abspath(self.get_output_directory())
