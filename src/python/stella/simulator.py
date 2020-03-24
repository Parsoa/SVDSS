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
    counttable,
    map_reduce,
    statistics,
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
        random.seed(c.seed)
        self.tracks = self.load_structural_variations()
        self.assign_event_zygosities()
        self.index_tracks()
        self.export_bed(self.tracks, 'all')
        self.export_bed(self.absent, 'absent')
        self.export_bed(self.present, 'present')
        self.export_bed(self.homozygous, 'homozygous')
        self.export_bed(self.heterozygous, 'heterozygous')
        self.extract_chromosomes()
        self.round_robin(self.chroms)

    def extract_chromosomes(self):
        self.chroms = extract_whole_genome()
        #for chrom, seq in self.chroms.items():
        #    print(chrom, 'with', len(seq), 'bases')
        #    with open(os.path.join(self.get_current_job_directory(), chrom + '.fa'), 'w') as chrom_file:
        #        chrom_file.write('>' + chrom + '\n')
        #        chrom_file.write(seq.strip())
        #        chrom_file.write('\n')

    def load_structural_variations(self):
        c = config.Configuration()
        tracks = [c.tracks[track] for track in c.tracks]
        return self.filter_overlapping_tracks(\
                    sorted(tracks, key = lambda x: (x.chrom, x.begin, x.end))
                )

    def filter_overlapping_tracks(self, tracks):
        remove = []
        i = 0
        while i < len(tracks):
            for j in range(i + 1, len(tracks)):
                # j is contained inside i
                if tracks[j].chrom != tracks[i].chrom:
                    i = j
                    break
                if tracks[j].begin <= tracks[i].end:
                    remove.append(j)
                    print(red(str(tracks[j])), 'overlaps', blue(str(tracks[i])))
                    continue
                else:
                    i = j
                    break
            if i == len(tracks) - 1:
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
        np.random.seed(92106429)
        for track in self.tracks:
            # psudo-deterministic STR on both strands
            random.seed(track.begin)
            # adjust sTR length and select repeat unit
            # TODO: maybe draw from a different ditribution
            #repeat_length = random.randint(2, 6)
            #track['repeat_length'] = repeat_length
            # generate repeat unit
            #repeat_unit = ''
            #bases = ['A', 'T', 'C', 'G']
            #for i in range(repeat_length):
            #    m = random.randint(0, 3)
            #    repeat_unit += bases[m]
            #track['repeat_unit'] = repeat_unit
            # generate alleles
            r = track.repeat
            num_repeats = int(r[:r.find('x')])
            repeat_unit = r[r.find('x') + 1:]
            track['repeat_unit'] = repeat_unit
            track['repeat_length'] = len(repeat_unit)
            p = np.random.poisson(num_repeats, 100)
            # reset random for allele selection
            random.seed(time.time())
            alleles = random.sample(list(p), k = 2)
            track['allele_father'] = alleles[0]
            track['allele_mother'] = alleles[1]
            track['seq'] = '.'
            # Make everything heterozygous for the sake of simplicity for now
            self.present.append(track)
        print(len(self.tracks), 'non-overlapping tracks')
        print(len(self.absent), '0|0')
        print(len(self.present), '1|1 or 1|0')
        print(len(self.homozygous), '1|1')
        print(len(self.heterozygous), '1|0')

    def index_tracks(self):
        self.track_index = {}
        for track in self.present:
            chrom = track.chrom.lower()
            if chrom not in self.track_index:
                self.track_index[chrom] = []
            self.track_index[chrom].append(track)
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
            if self.name == 'child':
                strand_1 = self.apply_events_to_chromosome('chry', 'father')
                strand_2 = self.apply_events_to_chromosome('chrx', 'mother')
        if chrom == 'chry':
            exit()
        self.export_diploid_chromosome_fasta(chrom, [strand_1, strand_2])
        self.export_fastq(seq, chrom + '_diploid')
        exit()

    def export_diploid_chromosome_fasta(self, chrom, strands):
        print('Exporting FASTA file:', green(chrom + '.fa')) 
        c = config.Configuration()
        with open(os.path.join(self.get_current_job_directory(), chrom + '_diploid.fa'), 'w') as fasta_file:
            for i in range(0, 2):
                seq = strands[i]
                fasta_file.write('>' + chrom + '_' + str(i + 1) + '\n')
                fasta_file.write(seq.strip())
                fasta_file.write('\n')
                continue
                l = len(seq)
                n = 100
                num_lines = l / n
                for i in range(0, num_lines):
                    line = seq[i * n : (i + 1) * n].upper() + '\n'
                    fasta_file.write(line)

    def export_chromosome_fasta(self, chrom, seq, name):
        print('Exporting FASTA file:', green(name + '.fa')) 
        c = config.Configuration()
        with open(os.path.join(self.get_current_job_directory(), name + '.fa'), 'w') as fasta_file:
            fasta_file.write('>' + chrom + '\n')
            fasta_file.write(seq.strip())
            fasta_file.write('\n')
            return
            l = len(seq)
            n = 100
            num_lines = l / n
            for i in range(0, num_lines):
                line = seq[i * n : (i + 1) * n].upper() + '\n'
                fasta_file.write(line)

    def apply_events_to_chromosome(self, chrom, haplotype):
        seq = ''
        previous = 0
        for track in self.track_index[chrom.lower()]:
            #if track.chrom.lower() != chrom.lower():
            #    print(yellow('skipping', track, 'on', chrom))
            #    continue
            #print('applying', track, 'to', chrom)
            seq += self.chroms[chrom][previous:track.begin]
            seq += track.repeat_unit * int(track['allele_' + haplotype])
            previous = track.end
        seq += self.chroms[chrom][previous:]
        return seq

    def export_fastq(self, seq, name):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        num_reads = len(seq) * c.coverage / 100
        fasta = os.path.join(self.get_current_job_directory(), name + '.fa')
        #command = "{} --data-type CCS --model_qc --accuracy-min 0.95 --length-max 15000 --length-min 5000  --length-sd 2000 --length-mean 10000 {}".format(pbsim, fasta)
        print(os.getcwd())
        command = "pbsim --data-type CCS --depth 30 --sample-fastq /share/hormozdiarilab/Codes/Fulminata/data/Samples/HG00731/HG00731_20190925_EEE_m54329U_190528_231241.Q20.fastq --prefix {} --length-max 20000 --length-min 2000 {}".format(name, fasta)
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT, cwd = self.get_current_job_directory())

    def reduce(self):
        c = config.Configuration()
        #self.export_reference_genome()
        #self.export_reference_jellyfish_table()
        return self.present

    def export_reference_genome(self):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        command = 'cat '
        for chrom in self.chroms:
            path = os.path.join(self.get_current_job_directory(), chrom + '.fa')
            command += path + ' ' 
        path = os.path.join(self.get_current_job_directory(), 'reference.fa')
        command += '> ' + path
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def export_reference_jellyfish_table(self):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        print('Generating Jellyfish table')
        command = "jellyfish count -m " + str(c.ksize) + " -s 1000000000 -t 24 --canonical --out-counter-len 2 "
        command += os.path.join(self.get_current_job_directory(), 'reference.fa')
        command += ' -o ' + os.path.join(self.get_current_job_directory(), 'reference_' + str(c.ksize) + 'k.jf')
        print(command)
        output = subprocess.call(command, shell = True, stdout = FNULL, stderr = subprocess.STDOUT)

    def export_bam(name, ref):
        c = config.Configuration()
        FNULL = open(os.devnull, 'w')
        fasta = os.path.join(get_output_directory(), 'reference.fa')
        fastq = os.path.join(get_output_directory(), 'test.fq')
        # generate sam file
        print('Generating SAM file')
        sam_file = os.path.join(get_output_directory(), 'test.sam')
        command = "bwa mem -M -t 24 {} {} > {}".format(fasta, fastq, sam_file)
        subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        # generate bam file
        print('Generating unsorted BAM file')
        unsorted_bam = os.path.join(get_output_directory(), name + '.unsorted.bam')
        command = "samtools view -S -b {} > {}".format(sam_file, unsorted_bam)  
        subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        #print('Sorting ...')
        #bam = os.path.join(get_output_directory(), name)
        #command = "samtools sort {} -o {}".format(unsorted_bam, bam)
        #subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        #print('Indexing ...')
        #bam_index = os.path.join(get_output_directory(), name + '.bai')
        #command = "samtools index {} {}".format(bam + '.bam', bam_index)
        #subprocess.call(command, shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
        print('Done!')

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
        self.index_tracks()
        self.extract_chromosomes()
        self.round_robin(self.chroms)

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
        #extract_whole_genome()
        self.tracks = {'father': True, 'mother': True}
        self.round_robin(self.tracks)

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
        np.random.seed(92106429)
        for parent in ['father', 'mother']:
            path = os.path.join(self.get_current_job_directory(), parent, 'present.bed')
            _tracks[parent] = bed.load_tracks_from_file(path)
        for f, m in zip(_tracks['father'], _tracks['mother']):
            t = copy.deepcopy(f)
            a = np.random.poisson(0, 4) - 5
            # assume the child is always male and naturally inherits father's allele
            # chrX is always inherited from mother
            if f.chrom.lower() == 'chrx':
                t['allele_father'] = '.'
                t['allele_mother'] = random.choice([int(m['allele_father']) + a[2], int(m['allele_mother']) + a[3]])
            # chrY always inherited from father and has a single haplotype
            elif f.chrom.lower() == 'chry':
                t['allele_father'] = int(f['allele_father']) + a[0]
                t['allele_mother'] = '.'
            # other chromosomes are normal
            else:
                t['allele_father'] = random.choice([int(f['allele_father']) + a[0], int(f['allele_mother']) + a[1]])
                t['allele_mother'] = random.choice([int(m['allele_father']) + a[2], int(m['allele_mother']) + a[3]])
            tracks.append(t)
        job = ChildSimulation(name = 'child', present = tracks)
        tracks = job.execute()

