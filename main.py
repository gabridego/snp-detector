import argparse
import os
import pickle
import time
from typing import Dict

import Bio.SeqIO

from utils import plot_frequency


def parse_args():
    parser = argparse.ArgumentParser(description="SNP detector")
    parser.add_argument('-p', '--path', nargs=2, help='Paths to FASTA files (or to stored binary files if -l)')
    parser.add_argument('-k', type=int, nargs='?', default=20, help='Length of patterns')
    parser.add_argument('-t', '--threshold', type=int, nargs='?', default=1, help='Threshold for K-mers filtering')
    parser.add_argument('-v', '--visualize', default=False, action='store_true', help='Plot intermediate results')
    parser.add_argument('-s', '--save', default=False, action='store_true', help='Save collected kmers')
    parser.add_argument('-l', '--load', default=False, action='store_true', help='Load collected kmers')
    return parser.parse_args()


def verify_args(args):
    if args.save and args.load:
        raise ValueError("-s and -l cannot be used at the same time")

    if args.load and (args.path[0].endswith('.fna') or args.path[1].endswith('.fna')):
        raise ValueError("in load mode, program expects binary files, not FASTA")


def parse_kmers(filename: str, k: int):
    kmers = {}

    print(f"start reading {os.path.basename(filename)}...")
    start = time.time()
    for record in Bio.SeqIO.parse(filename,
                                  "fasta"):
        seq = str(record.seq)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if kmer not in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
    end = time.time()

    print(f"{len(kmers)} {k}-mers extracted in {round(end - start, 2)} seconds.")

    return kmers


def filter_kmers(kmers: Dict[str, int], threshold: int = 1):
    filtered_kmers = {k: v for k, v in kmers.items() if v > threshold}
    kmers.clear()
    kmers.update(filtered_kmers)
    print(f"{len(kmers)} {len(list(kmers.keys())[0])}-mers kept after filtering (threshold={threshold})")


def main(args):
    """
    Collect patterns of length K, counting them in an hash table (dictionary)
    """
    k = args.k

    # read and filter wild strain
    if args.load:
        print(f"Load wild kmers from {args.path[0]}...")
        with open(args.path[0], 'rb') as f:
            wild_kmers = pickle.load(f)
    else:
        wild_kmers = parse_kmers(args.path[0], k)
        if args.save:
            filename = os.path.splitext(args.path[0])[0] + '_' + str(args.k) + '.pickle'
            print(f"Store wild kmers in {filename}...")
            with open(filename, 'wb') as f:
                pickle.dump(wild_kmers, f)

    if args.visualize:
        plot_frequency(wild_kmers, "Distribution of K-mers' number of occurrences for wild strain")

    filter_kmers(wild_kmers, args.threshold)

    if args.visualize:
        plot_frequency(wild_kmers, "Distribution of K-mers' number of occurrences for wild strain,"
                                   "after error filtering")

    # read and filter mutated strain
    if args.load:
        print(f"Load mutated kmers from {args.path[1]}...")
        with open(args.path[1], 'rb') as f:
            mut_kmers = pickle.load(f)
    else:
        mut_kmers = parse_kmers(args.path[1], k)
        if args.save:
            filename = os.path.splitext(args.path[1])[0] + '_' + str(args.k) + '.pickle'
            print(f"Store mutated kmers in {filename}...")
            with open(filename, 'wb') as f:
                pickle.dump(mut_kmers, f)

    if args.visualize:
        plot_frequency(mut_kmers, "Distribution of K-mers' number of occurrences for mutated strain")

    filter_kmers(mut_kmers, args.threshold)

    if args.visualize:
        plot_frequency(mut_kmers, "Distribution of K-mers' number of occurrences for mutated strain,"
                                  "after error filtering")


if __name__ == '__main__':
    args = parse_args()
    verify_args(args)

    main(args)
