import argparse
import os
import pickle
import time
from typing import Dict, List

import Bio.SeqIO
from colorama import init
from termcolor import colored

from utils import plot_frequency, levenshtein_distance


def parse_args():
    parser = argparse.ArgumentParser(description="SNP detector")
    parser.add_argument('-p', '--path', nargs=2,
                        help='Paths to FASTA files (or to stored binary files if -l)')
    parser.add_argument('-k', type=int, nargs='?', default=20,
                        help='Length of patterns')
    parser.add_argument('-f', '--filtering-threshold', type=int, nargs='?', default=1,
                        help='Threshold for K-mers filtering')
    parser.add_argument('-d', '--distance-threshold', type=int, nargs='?', default=10,
                        help='Threshold for Levenshtein distance')
    parser.add_argument('-v', '--visualize', default=False, action='store_true',
                        help='Plot intermediate results')
    parser.add_argument('-s', '--save', default=False, action='store_true',
                        help='Save collected kmers')
    parser.add_argument('-l', '--load', default=False, action='store_true',
                        help='Load collected kmers')
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


def get_unique_kmers(kmers1: Dict[str, int], kmers2: Dict[str, int], k: int, sort: bool = True):
    # build reads from kmers that are only in one of the two strains
    # this should allow us to keep some sequencing errors, at the end kmers with error will not be merged and can be
    # deleted safely
    # only check first and last k-1 characters, to deal with already merged sequences
    only1 = [kmer for kmer in kmers1 if kmer not in kmers2]
    only2 = [kmer for kmer in kmers2 if kmer not in kmers1]

    for only in [only1, only2]:
        z = ""
        # z is None when no sequence can be merged
        while z is not None:
            z = None
            for x in only:
                for y in only:
                    if x != y:
                        # y stars with x
                        if x[:k - 1] == y[len(y) - k + 1:]:
                            z = y + x[k - 1:]
                        # x starts with y
                        elif x[len(x) - k + 1:] == y[:k - 1]:
                            z = x + y[k - 1:]

                        # if the two sequences are concatenated, they can be removed from the list and the resulting
                        # sequence is inserted
                        if z is not None:
                            only.remove(x)
                            only.remove(y)
                            # insert at the beginning, may be more efficient
                            only.insert(0, z)
                            break
                if z is not None:
                    break

    if sort:
        only1.sort(reverse=True, key=lambda x: len(x))
        only2.sort(reverse=True, key=lambda x: len(x))

    return only1, only2


def print_snps(sequences1: List[str], sequences2: List[str], threshold: int = 10):
    print()
    print("Detected SNPs:")
    for x in sequences1:
        best_dist = float('inf')
        for y in sequences2:
            dist = levenshtein_distance(x, y)
            if dist < best_dist:
                best_string = y
                best_dist = dist

        if best_dist < threshold:
            snp_pos = []
            for i, (c1, c2) in enumerate(zip(x, best_string)):
                if c1 != c2:
                    snp_pos.append(i)

            for i, c in enumerate(x):
                if i in snp_pos:
                    print(colored(c, 'red'), end='')
                else:
                    print(c, end='')
            print(" -->")
            for i, c in enumerate(best_string):
                if i in snp_pos:
                    print(colored(c, 'red'), end='')
                else:
                    print(c, end='')
            print(f"\ndistance = {best_dist}")
            print()


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

    filter_kmers(wild_kmers, args.filtering_threshold)

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

    filter_kmers(mut_kmers, args.filtering_threshold)

    if args.visualize:
        plot_frequency(mut_kmers, "Distribution of K-mers' number of occurrences for mutated strain,"
                                  "after error filtering")

    only_wild, only_mut = get_unique_kmers(wild_kmers, mut_kmers, k)

    print_snps(only_wild, only_mut, args.distance_threshold)


if __name__ == '__main__':
    init()
    args = parse_args()
    verify_args(args)

    main(args)
