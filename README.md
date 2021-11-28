# Salmonella Outbreak Project

## Instructions

Detected K-mers are stored in a hash table (python dictionary) for each strain, where the key is the k-mer and the
value is the number of occurrences in the full file. The number of occurrences of each k-mer can therefore be accessed
in constant time.

``
python main.py -p data/salmonella-enterica.reads.fna data/salmonella-enterica-variant.reads.fna -k 50 -t 10 -v
``

### Command line interface

```
usage: main.py [-h] [-p PATH PATH] [-f [FORMAT]] [-k [K]]
               [-t [FILTERING_THRESHOLD]] [-d [DISTANCE_THRESHOLD]] [-v] [-s]
               [-l]

SNP detector

optional arguments:
  -h, --help            show this help message and exit
  -p PATH PATH, --path PATH PATH
                        Paths to FASTA files (or to stored binary files if -l)
  -f [FORMAT], --format [FORMAT]
                        Sequencing file format for Biopython
  -k [K]                Length of k-mers
  -t [FILTERING_THRESHOLD], --filtering-threshold [FILTERING_THRESHOLD]
                        Threshold for k-mers filtering
  -d [DISTANCE_THRESHOLD], --distance-threshold [DISTANCE_THRESHOLD]
                        Threshold for Levenshtein distance
  -v, --visualize       Plot intermediate results
  -s, --save            Save collected k-mers
  -l, --load            Load collected k-mers
```

If the `filtering-threshold` argument is not provided, user is interactively asked to input a value during execution. 

### Storing and loading dictionaries
Interface is provided to store/load the dictionary of detected k-mers in/from binary files using
[pickle](https://docs.python.org/3/library/pickle.html). This allows to test different thresholds for the filters,
detecting the k-mers only once and saving time.

**Warning:** the resulting binary files can be huge.

#### Store
To store the dictionaries computed in the run:

``
python main.py -p data/salmonella-enterica.reads.fna data/salmonella-enterica-variant.reads.fna -k 20 -v -s
``

This will create two binary files, called `data/salmonella-enterica.reads_20.pickle` and
`data/salmonella-enterica-variant.reads_20.pickle` (`20` as the provided `k`).

#### Load
To load a previously stored binary file:

``
python main.py -p data/salmonella-enterica.reads_20.pickle data/salmonella-enterica-variant.reads_20.pickle -k 20 -v -l
``

## Data

Data is not included in this repo, please download it from the
[course website](https://clovisg.github.io/teaching/protein-structure-prediction/sequences/). A sample file for
testing can be found
[here](https://gitlab.ensimag.fr/galiezc/protein-structure-prediction/-/blob/master/hands-on/reference-data/salmonella-enterica.reads.fna).

COVID-19 sequences can be downloaded from the [COVID-19 Data Portal](https://www.covid19dataportal.org/sequences),
looking for entries for which raw reads are available. For example, Illumina reads for lineages B.1.1.7 and B.1.1.8 can
be respectively downloaded [here](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR468/008/ERR4682028/ERR4682028_1.fastq.gz) and
[here](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR445/003/ERR4458283/ERR4458283_1.fastq.gz).

## Authors
The project is presented by the CEO, CTO and CHO of DZA Computing:
- Sophie Zhang
- Enrico Agrippino
- Gabriele Degola

© 2021 DZA Computing. All rights reserved.