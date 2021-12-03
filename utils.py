from typing import Dict, Optional

from matplotlib import pyplot as plt
import numpy as np


def plot_frequency(kmers: Dict[str, int], title: Optional[str] = None, log: bool = False):
    counts = list(kmers.values())
    _, ax = plt.subplots(figsize=(20, 15))
    ax.hist(counts, bins=len(set(counts)), log=log)
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.set_xlabel("Number of occurrences")
    ax.set_ylabel("Frequency")
    if title:
        ax.set_title(title)
    plt.show()


def levenshtein_distance(token1: str, token2: str):
    # https://blog.paperspace.com/implementing-levenshtein-distance-word-autocomplete-autocorrect/
    distances = np.zeros((len(token1) + 1, len(token2) + 1))

    for t1 in range(len(token1) + 1):
        distances[t1][0] = t1

    for t2 in range(len(token2) + 1):
        distances[0][t2] = t2
    
    for t1 in range(1, len(token1) + 1):
        for t2 in range(1, len(token2) + 1):
            if token1[t1 - 1] == token2[t2 - 1]:
                distances[t1][t2] = distances[t1 - 1][t2 - 1]
            else:
                a = distances[t1][t2 - 1]
                b = distances[t1 - 1][t2]
                c = distances[t1 - 1][t2 - 1]
                
                if a <= b and a <= c:
                    distances[t1][t2] = a + 1
                elif b <= a and b <= c:
                    distances[t1][t2] = b + 1
                else:
                    distances[t1][t2] = c + 1

    return distances[len(token1)][len(token2)]