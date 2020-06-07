# SeqMutations Utils

### Description
Utility to generate codon mutations of a given coding DNA sequence.


### Usage
```text
usage: seq_mutations_utils.py [-h] -i INPUT_FILE -p INITIAL_POSITION
                              [-m MUTATIONS] -n NUM -o OUTPUT_PREFIX [-a] [-v]
                              [-x]

Mutations in substring

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file INPUT_FILE
                        File containing the coding DNA sequence
  -p INITIAL_POSITION, --initial-pos INITIAL_POSITION
                        Initial position of substring to modify starting from
                        1
  -m MUTATIONS, --mutations MUTATIONS
                        List of mutations in the form of
                        '(ala,glu,...),(val,...),...' -- NOTE: the "'" must be
                        included
  -n NUM, --num NUM     Number of consecutive codons
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output prefix file name
  -a, --all-codons      Try all codons
  -v, --inv             Also consider combinations with original codons
                        without mutations
  -x, --single-output-file
                        Create single fasta file
```
