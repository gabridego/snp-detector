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
filt1=filter_kmers(kmers1,threshold=18)
filt2=filter_kmers(kmers2,threshold=18)

plot_frequency(filt1, "Distribution of K-mers' number of occurrences for wild strain,"
                           "after error filtering")

import matplotlib.pyplot as plt
counts = list(kmers1.values())
fig, ax = plt.subplots(figsize=(20, 15))
ax.hist(counts, bins=len(set(counts)))
ax.set_xlabel("Number of occurrences")
ax.set_ylabel("Frequency")
ax.set_xlim([0, 25])
ax.set_ylim([0, 25000])
plt.show()

plot_frequency(filt2, "Distribution of K-mers' number of occurrences for mutated strain,"
                           "after error filtering")

sum(x > 100 for x in kmers1.values())
sum(x > 1000 for x in kmers1.values())
sum(x < 15 and x>10 for x in kmers1.values())

sum(x not in filt2 for x in filt1.keys())


