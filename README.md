# Salmonella Outbreak Project

## Instructions

With `hash` option, store patterns in an hash table (python dictionary), where the key is the sequence of nucleotides
and the value is the number of occurrences in the full file. The value of each pattern can be accessed in constant time.

With `tree` option, store patterns in four trees, where the roots are the existing nucleotides. The depth of the trees
is always `K`. This is computationally expensive.

``
python main.py -p salmonella-enterica.reads.fna -t <hash/tree> -k 30
``

## Authors
Enrico Agrippino  
Gabriele Degola  
Sophie Zhang
