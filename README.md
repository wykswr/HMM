# HMM Viterbi Sequence Aligner

## Algorithm
This model utilizes paired hidden markov model (HMM) to find the alignment between 2 sequences.

## Usage
The model is implemented as a UNIX command line tool, with syntax:

`HMM -i <infile> -a <path output> -b <alignment output>`

in which "path output" contains the hidden states path, "alignment output" contains alignment results.

## Compilation
cmake minimum required VERSION: 3.23