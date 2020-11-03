import io
import os
import pwd
import sys
import json
import time
import subprocess

from stella import (
    config,
)

from stella.logger import *

print = pretty_print

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

class BedTrack:

    def __init__(self, chrom, begin, end, fields = []):
        self.chrom = chrom
        self.begin = begin
        self.end = end
        self.fields = []
        for (k, v) in fields:
            self[k] = v
        if not hasattr(self, 'svtype'):
            self.svtype = 'DEL'
        if hasattr(self, 'genotype'):
            if self.genotype == '0/1' or self.genotype == './1' or self.genotype == '1/.':
                self.genotype = '1/0'
        if not hasattr(self, 'id'):
            self.id = self.svtype + '@' + self.chrom + '_' + str(self.begin) + '_' + str(self.end)

    def lift(self):
        with open('liftover_tmp.bed', 'w') as b:
            b.write(self.chrom + '\t' + str(self.begin) + '\t' + str(self.end) + '\n')
        command = '/home/pkhorsand/local/bin/liftOver ' + 'liftover_tmp.bed' + ' /afs/genomecenter.ucdavis.edu/home/pkhorsand/hg19ToHg38.over.chain liftover_res.bed liftover_un.bed'
        output = subprocess.call(command, shell = True)
        track = load_tracks_from_file('liftover_res.bed')[0]
        os.remove('liftover_tmp.bed')
        os.remove('liftover_res.bed')
        os.remove('liftover_un.bed')
        return track

    @staticmethod
    def json_serialize(self):
        d = {}
        d['end'] = self.end
        d['begin'] = self.begin
        d['chrom'] = self.chrom
        for (k, v) in self.fields:
            d[k] = v
        return d

    def serialize(self):
        s = ''
        s += str(self.chrom) + '\t' + str(self.begin) + '\t' + str(self.end)
        for key in self.fields:
            s += '\t' + str(self[key])
        s += '\n'
        return s

    def header(self):
        s = '#CHROM\tBEGIN\tEND'
        for key in self.fields:
            s += '\t' + key.upper()
        s += '\n'
        return s

    def __str__(self):
        return self.svtype + '@' + self.chrom + '_' + str(self.begin) + '_' + str(self.end)
        return self.id

    def __getitem__(self, key):
        return getattr(self, key.lower())

    def __setitem__(self, key, value):
        key = key.lower()
        setattr(self, key, value)
        if not key in self.fields:
            self.fields.append(key)

    def __delattr__(self, key):
        key = key.lower()
        del self.__dict__[key]
        self.fields.remove(key)

# ============================================================================================================================ #
# BED Tracks
# ============================================================================================================================ #

def as_dict(tracks):
    return {str(track): track for track in tracks}

# sorts a dictionary of tracks into a list
def sort_tracks(tracks):
    if isinstance(tracks, dict):
        _tracks = [tracks[track] for track in tracks]
    else:
        _tracks = tracks
    return sorted(_tracks, key = lambda x: (x.chrom, x.begin, x.end))
    #return sorted(sorted(_tracks, key = lambda x: x.begin), key = lambda y: y.chrom)

def filter_short_tracks(tracks):
    return [track for track in tracks if abs(int(track.svlen)) > 50]

def filter_overlapping_tracks(tracks, svtype):
    i = 0
    remove = []
    tracks = sort_tracks(tracks)
    while i < len(tracks):
        for j in range(i + 1, len(tracks)):
            if tracks[j].chrom != tracks[i].chrom:
                i = j
                break
            if svtype == 'DEL':
                if tracks[j].begin <= tracks[i].end:
                    remove.append(j)
                    user_print_warning(str(tracks[j]), white('overlaps'), blue(str(tracks[i])))
                    continue
            if svtype == 'INS':
                if tracks[j].begin - tracks[i].end < 100:
                    remove.append(j)
                    user_print_warning(str(tracks[j]), white('is too close to'), blue(str(tracks[i])))
                    continue
            if svtype == 'ALL':
                if tracks[j].begin <= tracks[i].end:
                    remove.append(j)
                    user_print_warning(str(tracks[j]), white('overlaps previous track'), blue(str(tracks[i])))
                    continue
                if tracks[j].begin - tracks[i].end < 500:
                    remove.append(j)
                    user_print_warning(str(tracks[j]), white('too close to previous track'), blue(str(tracks[i])))
                    continue
            i = j
        if j == len(tracks) - 1:
            break
    n = 0
    for index in sorted(remove):
        tracks.pop(index - n)
        n = n + 1
    return tracks

def track_from_id(name):
    svtype, coords = name.split('@')
    tokens = coords.split('_')
    return BedTrack(tokens[0], int(tokens[1]), int(tokens[2]), [('SVTYPE', svtype)])

# keywords is an array of tuples: (name, default, transformation)
def load_tracks_from_file(path, parse_header = True, keywords = []):
    tracks = []
    with open(path, 'r') as f:
        if parse_header:
            header = f.readline()
            fields = header.lower().split()
        line = f.readline()
        while line:
            tokens = line.split()
            kwargs = []
            for i in range(3, len(tokens)):
                kwargs.append((fields[i], tokens[i]))
            track = BedTrack(chrom = tokens[0], begin = int(tokens[1]), end = int(tokens[2]), fields = kwargs)
            tracks.append(track)
            line = f.readline()
    return tracks

def load_tracks_from_file_as_dict(path, parse_header = True, keywords = []):
    return {str(track): track for track in load_tracks_from_file(path, parse_header, keywords)}


