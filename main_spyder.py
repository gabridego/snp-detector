# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 17:49:02 2021

@author: enrya
"""
import argparse
import os
import time
from typing import Dict
import pickle

import Bio.SeqIO

from utils import plot_frequency

k=50
# with k=25 36197708 kmers1 and 36197979 kmers2
filename1="data/salmonella-enterica.reads.fna"
kmers1 = {}

with open("../data/salmonella-enterica_50.pickle", 'rb') as f:
    kmers1 = pickle.load(f)

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

with open("../data/salmonella-enterica-variant_50.pickle", 'rb') as f:
    kmers2 = pickle.load(f)

print(f"start reading {os.path.basename(filename2)}...")
start = time.time()
n_reads=0
for record in Bio.SeqIO.parse(filename2,
                              "fasta"):
    seq = str(record.seq)
    n_reads+=1
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer not in kmers2:
            kmers2[kmer] = 0
        kmers2[kmer] += 1
end = time.time()

with open("data/salmonella-enterica_50.pickle", 'wb') as f:
    pickle.dump(kmers1, f)

print(f"{len(kmers2)} {k}-mers extracted in {round(end - start, 2)} seconds.")

def filter_kmers(kmers: Dict[str, int], threshold=1):
    filtered_kmers = {k: v for k, v in kmers.items() if v > threshold}
    print(f"{len(filtered_kmers)} {len(list(kmers.keys())[0])}-mers kept after filtering (threshold={threshold})")
    return filtered_kmers

# plot entire strains
plot_frequency(kmers1, f"Distribution of K-mers' number of occurrences for wild strain,k={k}")
plot_frequency(kmers2, f"Distribution of K-mers' number of occurrences for mutate strain,k={k}")

#filtering
#filt1=filter_kmers(kmers1,threshold=8)
#filt2=filter_kmers(kmers2,threshold=14)
filt1=filter_kmers(kmers1,threshold=10)
filt2=filter_kmers(kmers2,threshold=8)

plot_frequency(filt1, f"Distribution of K-mers' number of occurrences for wild strain,"
                           f"after error filtering, k={k}, filt=8")

import matplotlib.pyplot as plt
counts = list(kmers1.values())
fig, ax = plt.subplots(figsize=(20, 15))
ax.hist(counts, bins=len(set(counts)))
ax.set_xlabel("Number of occurrences")
ax.set_ylabel("Frequency")
ax.set_xlim([0, 25])
ax.set_ylim([0,100])
ax.set_title("wild strain zoom, fliter=8, k=20")
plt.show()

import matplotlib.pyplot as plt
counts = list(kmers2.values())
fig, ax = plt.subplots(figsize=(20, 15))
ax.hist(counts, bins=len(set(counts)))
ax.set_xlabel("Number of occurrences")
ax.set_ylabel("Frequency")
ax.set_xlim([0, 25])
ax.set_ylim([0, 100])
ax.set_title("mutate strain zoom, fliter=8, k=20")
plt.show()

plot_frequency(filt2, f"Distribution of K-mers' number of occurrences for mutated strain,"
                           f"after error filtering, filt=8,k={k}")

sum(x == 2 for x in kmers1.values())
sum(x == 3 for x in kmers1.values())
sum(x == 4 for x in kmers1.values())
sum(x == 5 for x in kmers1.values())
sum(x == 6 for x in kmers1.values())
sum(x > 1000 for x in kmers1.values())
sum(x < 12 and x>3 for x in kmers1.values())
sum(x not in filt2 for x in filt1.keys())
sum(x not in filt1 for x in filt2.keys())

# k=20, f=12 ==> 85 and 65
# k=20, f=10 ==> 44 and 158
# k=20, f1= 10, f2=14 ==> 44 and 138
# k=20, f1=9, f2=14 ==> 44 and 105
# k=20, f1=8, f2=13 ==> 47 and 88

def opt_f(kmers1,kmers2):
    minimum=1000000
    f1_opt=0
    f2_opt=0
    f2_dict={}
    for f1 in range(3,16):
        filt1=filter_kmers(kmers1,threshold=f1)
        for f2 in range(3,16):
            if f2 not in f2_dict:
                f2_dict[f2]=filter_kmers(kmers2,threshold=f2)
            filt2=f2_dict[f2]
            np1=sum(x not in filt2 for x in filt1.keys())
            np2=sum(x not in filt1 for x in filt2.keys())
            if abs(np1-np2)+(np1+np2)<minimum:
                minimum=abs(np1-np2)+(np1+np2)
                f1_opt=f1
                f2_opt=f2
    return(f1_opt,f2_opt,minimum)

f1_opt,f2_opt,minimum=opt_f(kmers1,kmers2)
                
"""
Results of minimization:
    minimum=136
    f1_opt=8
    f2_opt=14
"""

## build reads from kmers that are only in one of the two strains
## this should allow us to keep some sequencing errors, at the end kmers with error will not be merged and can be deleted safely
## only check first and last k characters, to deal with already merged sequences
only1 = [kmer for kmer in filt1 if kmer not in filt2]
only2 = [kmer for kmer in filt2 if kmer not in filt1]

for only in [only1, only2]:
    z = ""
    # z is None when no more sequences can be merged
    while z is not None:
        z = None
        for x in only:
            for y in only:
                if x != y:
                    # y stars with x
                    if x[:k-1] == y[len(y)-k+1:]:
                        z = y + x[k-1:]
                    # x starts with y
                    elif x[len(x)-k+1:] == y[:k-1]:
                        z = x + y[k-1:]
                        
                    # if the two sequences are concatenated, they can be removed from the list and the resulting sequence is inserted
                    if z is not None:
                        only.remove(x)
                        only.remove(y)
                        # insert at the beginning, may be more efficient
                        only.insert(0, z)
                        break
            if z is not None:
                break

only1.sort(reverse=True, key=lambda x: len(x))
only2.sort(reverse=True, key=lambda x: len(x))

print(only1)
print(only2)

print("Detected SNPs:")
for x in only1:
    best_dist = float('inf')
    for y in only2:
        dist = levenshteinDistance(x, y)
        if dist < best_dist:
            best_string = y
            best_dist = dist
    
    if best_dist < 10:
        x = list(x)
        best_string = list(best_string)
        for i, (c1, c2) in enumerate(zip(x, best_string)):
            x[i] = c1
            best_string[i] = c2
            if c1 != c2:
                x[i] += "\u0332"
                best_string[i] += "\u0332"
        x = "".join(x)
        best_string = "".join(best_string)
        
        print(f"{x} -->\n{best_string}\ndistance = {best_dist}")
        print()
