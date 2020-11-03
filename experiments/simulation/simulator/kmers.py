from __future__ import print_function

import io
import os
import re
import pwd
import sys
import copy
import math
import time
import random

# ============================================================================================================================ #
# kmer helpers
# ============================================================================================================================ #

def canonicalize(seq):
    seq = seq.upper()
    reverse_complement = reverse_complement_sequence(seq)
    return seq if seq < reverse_complement else reverse_complement

def c_extract_kmers(k = 32, counter = lambda x: 1, count = 1, overlap = True, canonical = True, *args):
    kmers = {}
    for s in args:
        i = 0
        while i <= len(s) - k and i >= 0:
            kmer = canonicalize(s[i : i + k]) if canonical else s[i : i + k]
            cc = counter(kmer)
            if cc > count:
                i += 1
                continue
            if not kmer in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
            if not overlap:
                i += k
            else:
                i += 1
    return kmers

def extract_kmers(k = 32, canonical = True, *args):
    kmers = {}
    for s in args:
        for i in range(0, len(s) - k + 1):
            kmer = canonicalize(s[i : i + k]) if canonical else s[i : i + k]
            if not kmer in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
    return kmers

def stream_kmers(k = 32, canonical = True, overlap = True, *args):
    for s in args:
        for i in range(0, len(s) - k + 1, 1 if overlap else k):
            kmer = s[i : i + k]
            yield canonicalize(kmer) if canonical else kmer

def index_kmers(k = 32, canonical = True, *args):
    kmers = {}
    index = 0
    for kmer in stream_kmers(k, canonical, *args):
        if not kmer in kmers:
            kmers[kmer] = []
        kmers[kmer].append(index)
        index += 1
    return kmers

def find_kmer(k, kmers):
    rc = reverse_complement(k)
    if k in kmers:
        return k
    if rc in kmers:
        return rc
    return None

def reverse_complement(seq):
    return complement_sequence(seq[::-1])

def reverse_complement_sequence(seq):
    return complement_sequence(seq[::-1])

def complement_sequence(seq):
    # A-> C and C->A
    seq = seq.replace('A', 'Z')
    seq = seq.replace('T', 'A')
    seq = seq.replace('Z', 'T')
    #
    seq = seq.replace('G', 'Z')
    seq = seq.replace('C', 'G')
    seq = seq.replace('Z', 'C')
    #
    return seq

def is_subsequence(x, y):
    if len(y) < len(x):
        return False
    it = iter(y)
    #debug_log('Checking', green(x), 'substring of', blue(y))
    return all(c in it for c in x)

def is_canonical_subsequence(x, y):
    if len(y) < len(x):
        return False
    #debug_log('Checking', green(x), 'substring of', blue(y))
    return is_subsequence(reverse_complement(x), y) or is_subsequence(x, y)
    #return all(c in it for c in x) or all(c in it for c in reverse_complement(x))

def find_all(string, substring):
    l = []
    index = -1
    while True:
        index = string.find(substring, index + 1)
        if index == -1:
            break
        l.append(index)
    return l

def calculate_gc_content(seq):
    n = 0
    for c in seq.upper():
        if c == 'C' or c == 'G':
            n += 1
    return n
    #return len(list(filter(lambda x: x == 'G' or x == 'C', seq)))

def is_kmer_low_entropy(kmer):
    kmers = {}
    m = 0
    for k in stream_kmers(3, False, True, kmer):
        if not k in kmers:
            kmers[k] = 0
        kmers[k] += 1
        if kmers[k] > m:
            m = kmers[k]
    if len(kmers) < 7:
        return True
    if m > (len(kmer) - 3 + 1) / 3:
        return True
    return False
