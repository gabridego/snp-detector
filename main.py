import argparse
import time

import Bio.SeqIO


class Node(object):
    """
    Node object for tree building
    """
    def __init__(self, data):
        self.data = data
        self.children = {}
        self.parent = None
        self.count = 0

    def __str__(self):
        return f"{self.data}: parent {self.parent.data}, children {list(self.children.keys())}"

    def add_child(self, obj):
        self.children[obj.data] = obj

    def set_parent(self, obj):
        self.parent = obj

    def get_child(self, data):
        return self.children[data]


def parse_args():
    parser = argparse.ArgumentParser(description="Salmonella outbreak")
    parser.add_argument('-p', '--path', nargs='+', help='Paths to FASTA files (max 2)')
    parser.add_argument('-t', '--type', nargs='?', choices=['tree', 'hash'],
                        default='hash', help='What to do with the files')
    parser.add_argument('-k', type=int, nargs='?', default=30, help='Length of patterns')
    return parser.parse_args()


def use_trees(args):
    """
    Build tree for nucleotide patterns, of depth K
    """
    seqs = []

    for record in Bio.SeqIO.parse(args.path[0],
                                  "fasta"):
        seqs.append(str(record.seq))

    roots = {}
    K = args.k
    leaves = set()

    start = time.time()
    for seq in seqs:
        for i in range(len(seq) - K + 1):
            s = seq[i:i + K]
            if s[0] not in roots:
                roots[s[0]] = Node(s[0])
            parent = roots[s[0]]
            for j in range(1, len(s)):
                if s[j] not in parent.children:
                    new_node = Node(s[j])
                    parent.add_child(new_node)
                    new_node.set_parent(parent)
                parent = parent.get_child(s[j])
            parent.count += 1
            leaves.add(parent)
    end = time.time()

    print(f"tree built in {round(end - start, 2)} seconds.")

    return roots


def use_hash(args):
    """
    Collect patterns of length K, counting them in an hash table (dictionary)
    """
    hash_tab = {}
    K = args.k

    start = time.time()
    for record in Bio.SeqIO.parse(args.path[0],
                                  "fasta"):
        seq = str(record.seq)
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i + K]
            if kmer not in hash_tab:
                hash_tab[kmer] = 0
            hash_tab[kmer] += 1
    end = time.time()

    print(f"{len(hash_tab)} patterns of {K} characters collected in {round(end - start, 2)} seconds.")


if __name__ == '__main__':
    args = parse_args()

    if args.type == 'hash':
        use_hash(args)
    else:
        use_trees(args)
