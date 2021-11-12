# Salmonella Outbreak Project

## Instructions

Detected K-mers are stored in an hash table (python dictionary) for each strain, where the key is the K-mer and the
value is the number of occurrences in the full file. The number of occurrences of each K-mer can therefore be accessed
in constant time.

``
python main.py -p data/salmonella-enterica.reads.fna data/salmonella-enterica-variant.reads.fna -k 20 -f 5 -v
``

### Command line interface

```
usage: main.py [-h] [-p PATH PATH] [-k [K]] [-f [FILTERING_THRESHOLD]]
               [-d [DISTANCE_THRESHOLD]] [-v] [-s] [-l]

SNP detector

optional arguments:
  -h, --help            show this help message and exit
  -p PATH PATH, --path PATH PATH
                        Paths to FASTA files (or to stored binary files if -l)
  -k [K]                Length of patterns
  -f [FILTERING_THRESHOLD], --filtering-threshold [FILTERING_THRESHOLD]
                        Threshold for K-mers filtering
  -d [DISTANCE_THRESHOLD], --distance-threshold [DISTANCE_THRESHOLD]
                        Threshold for Levenshtein distance
  -v, --visualize       Plot intermediate results
  -s, --save            Save collected kmers
  -l, --load            Load collected kmers
```

### Storing and loading dictionaries
Interface is provided to store/load the dictionary of detected K-mers in/from binary files using
[pickle](https://docs.python.org/3/library/pickle.html). his allows to test different thresholds for the filters,
detecting the K-mers only once and saving time.

**Warning:** the resulting binary files can be huge.

#### Store
To store the dictionaries computed in the run:

``
python main.py -p data/salmonella-enterica.reads.fna data/salmonella-enterica-variant.reads.fna -k 20 -v -s
``

This will create two binary files, called `data/salmonella-enterica.reads_20.fna` and
`data/salmonella-enterica-variant.reads_20.fna` (`20` as the provided `k`).

#### Load
To load a previously stored binary file:

``
python main.py -p data/salmonella-enterica_20.reads.fna data/salmonella-enterica-variant_20.reads.fna -k 20 -v -l
``

## Data

Data is not included in this repo, please download it from the
[course website](https://clovisg.github.io/teaching/protein-structure-prediction/sequences/). A sample file for
testing can be found
[here](https://gitlab.ensimag.fr/galiezc/protein-structure-prediction/-/blob/master/hands-on/reference-data/salmonella-enterica.reads.fna).

## Authors
The project is presented by the CEO, CTO and CHO of DZA Computing:
- Sophie Zhang
- Enrico Agrippino
- Gabriele Degola

© 2021 DZA Computing. All rights reserved.