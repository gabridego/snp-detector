from typing import Dict

from matplotlib import pyplot as plt


def plot_frequency(kmers: Dict[str, int], title=None):
    counts = list(kmers.values())
    _, ax = plt.subplots(figsize=(20, 15))
    ax.hist(counts, bins=len(set(counts)))
    ax.set_xlabel("Number of occurrences")
    ax.set_ylabel("Frequency")
    if title:
        ax.set_title(title)
    plt.show()
