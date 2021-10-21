import argparse
import os
import time
from typing import Dict

import Bio.SeqIO

from utils import plot_frequency


def parse_args():
    parser = argparse.ArgumentParser(description="SNP detector")
    parser.add_argument('-p', '--path', nargs=2, help='Paths to FASTA files')
    parser.add_argument('-k', type=int, nargs='?', default=20, help='Length of patterns')
    parser.add_argument('-v', '--visualize', default=False, action='store_true', help='Plot intermediate results')
    return parser.parse_args()


def parse_kmers(filename: str, k: int):
    kmers = {}

    print(f"start reading {os.path.basename(args.path[0])}...")
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

    print(f"{len(kmers)} patterns of {k} characters extracted in {round(end - start, 2)} seconds.")

    return kmers


def filter_kmers(kmers: Dict[str, int], threshold=1):
    filtered_kmers = {k: v for k, v in kmers.items() if v > threshold}
    kmers.clear()
    kmers.update(filtered_kmers)
    print(f"{len(kmers)} patterns kept after filtering (threshold={threshold})")


def main(args):
    """
    Collect patterns of length K, counting them in an hash table (dictionary)
    """
    k = args.k

    # read and filter wild strain
    wild_kmers = parse_kmers(args.path[0], k)

    if args.visualize:
        plot_frequency(wild_kmers, "Distribution of K-mers' number of occurrences for wild strain")

    filter_kmers(wild_kmers)

    if args.visualize:
        plot_frequency(wild_kmers, "Distribution of K-mers' number of occurrences for wild strain,"
                                   "after error filtering")

    # read and filter mutated strain
    mut_kmers = parse_kmers(args.path[1], k)

    if args.visualize:
        plot_frequency(mut_kmers, "Distribution of K-mers' number of occurrences for mutated strain")

    filter_kmers(mut_kmers)

    if args.visualize:
        plot_frequency(mut_kmers, "Distribution of K-mers' number of occurrences for mutated strain,"
                                  "after error filtering")

    # compare the k-mers collected for the two strains


if __name__ == '__main__':
    args = parse_args()

    main(args)
