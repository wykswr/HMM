# HMM Viterbi Sequence Aligner

## Algorithm
This model utilizes paired hidden markov model (HMM) to find the alignment between 2 sequences.

## Usage
The model is implemented as a UNIX command line tool, with syntax:

`HMM [-f] -i <infile> -a <path | probability> -b <alignment>`

* `-f`: Add this flag to enable forward mode.
* `infile`: The input fastq file, in which every 2 sequences are considered as a pair.
* `path | probability`: The output path of inferred hidden states or probability in forward mode.
* `alignment`: The output path of alignment result, and is useless in forward mode.

## Compilation
cmake minimum required VERSION: 3.23