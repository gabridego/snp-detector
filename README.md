# Salmonella Outbreak Project

## Instructions

Detected K-mers are stored in an hash table (python dictionary) for each strain, where the key is the K-mer and the
value is the number of occurrences in the full file. The number of occurrences of each K-mer can therefore be accessed
in constant time.

``
python main.py -p data/salmonella-enterica.reads.fna data/salmonella-enterica-variant.reads.fna -k 20 -v
``

Data is not included in this repo, please download it from the
[course website](https://clovisg.github.io/teaching/protein-structure-prediction/sequences/).

## Authors
The project is presented by the CEO, CTO and CHO of DZA Computing:
- Sophie Zhang
- Enrico Agrippino
- Gabriele Degola

© 2021 DZA Computing. All rights reserved.