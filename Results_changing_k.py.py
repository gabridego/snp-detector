# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 15:54:20 2021

@author: enrya
"""
import numpy as np
from typing import Dict
import Bio.SeqIO

a=2 #starting value for f1 and f2
b=20 # stopping value for f1 and f2
filename1="Salmonella outbreak/data/salmonella-enterica.reads.fna"
filename2="Salmonella outbreak/data/salmonella-enterica-variant.reads.fna"

def filter_kmers(kmers: Dict[str, int], threshold=1):
    filtered_kmers = {k: v for k, v in kmers.items() if v > threshold}
    #print(f"{len(kmers)} {len(list(filtered_kmers.keys())[0])}-mers kept after filtering (threshold={threshold})")
    return filtered_kmers

def reading_kmers(filename):
    kmers = {}
    for record in Bio.SeqIO.parse(filename1,
                                  "fasta"):
        seq = str(record.seq)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if kmer not in kmers1:
                kmers1[kmer] = 0
            kmers1[kmer] += 1
    return kmers
    
for k in range(14,26):    # for on the reading of k-mers
    print(f"Reading {k}mers...")
    kmers1=reading_kmers(filename1)
    kmers2=reading_kmers(filename2)

    # Now for cycle on the filter thresholds
    f2_dict={}
    Not_present_1=np.zeros([b+1-a,b+1-a],dtype=int)
    Not_present_2=np.zeros([b+1-a,b+1-a],dtype=int)
    for f1 in range(a,b+1):
        filt1=filter_kmers(kmers1,threshold=f1)
        for f2 in range(a,b+1):
            if f2 not in f2_dict:
                f2_dict[f2]=filter_kmers(kmers2,threshold=f2)
            filt2=f2_dict[f2]
            Not_present_1[f1-a,f2-a]=sum(x not in filt2 for x in filt1.keys())
            Not_present_2[f1-a,f2-a]=sum(x not in filt1 for x in filt2.keys())
            
            # Saving results for the fixed k, f1 and f2 into a file
            print(f"Saving results for k={k} f1={f1} f2={f2}")
            np.savetxt(f"kmers_missings/np1_k{k}_f1_{f1}_f2_{f2}", Not_present_1)
            np.savetxt(f"kmers_missings/np2_k{k}_f1_{f1}_f2_{f2}", Not_present_2)
            
# trial on how to load file, without re running all the code
#prova=np.loadtxt(results_filename1,dtype=int)
