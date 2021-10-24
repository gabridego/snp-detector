# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 17:49:02 2021

@author: enrya
"""
import argparse
import os
import time
from typing import Dict

import Bio.SeqIO

from utils import plot_frequency

k=20
filename1="data/salmonella-enterica.reads.fna"
kmers1 = {}

print(f"start reading {os.path.basename(filename1)}...")
start = time.time()
for record in Bio.SeqIO.parse(filename1,
                              "fasta"):
    seq = str(record.seq)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer not in kmers1:
            kmers1[kmer] = 0
        kmers1[kmer] += 1
end = time.time()

print(f"{len(kmers1)} {k}-mers extracted in {round(end - start, 2)} seconds.")

filename2="data/salmonella-enterica-variant.reads.fna"
kmers2 = {}

print(f"start reading {os.path.basename(filename2)}...")
start = time.time()
for record in Bio.SeqIO.parse(filename2,
                              "fasta"):
    seq = str(record.seq)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer not in kmers2:
            kmers2[kmer] = 0
        kmers2[kmer] += 1
end = time.time()

print(f"{len(kmers2)} {k}-mers extracted in {round(end - start, 2)} seconds.")

def filter_kmers(kmers: Dict[str, int], threshold=1):
    filtered_kmers = {k: v for k, v in kmers.items() if v > threshold}
    print(f"{len(kmers)} {len(list(kmers.keys())[0])}-mers kept after filtering (threshold={threshold})")
    return filtered_kmers

# plot entire strains
plot_frequency(kmers1, "Distribution of K-mers' number of occurrences for wild strain")
plot_frequency(kmers2, "Distribution of K-mers' number of occurrences for mutate strain")

#filtering
filt1=filter_kmers(kmers1,threshold=5)
filt2=filter_kmers(kmers2,threshold=5)

plot_frequency(filt1, "Distribution of K-mers' number of occurrences for wild strain,"
                           "after error filtering")
plot_frequency(filt2, "Distribution of K-mers' number of occurrences for mutated strain,"
                           "after error filtering")

import numpy as np
#np.sum(np.array(kmers1.values()) > 100)
