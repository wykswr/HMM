# HMM Viterbi Sequence Aligner

## Algorithm

This model utilizes paired hidden markov model (HMM) to find the alignment between 2 sequences.

## Usage

The model is implemented as a UNIX command line tool, with syntax:

`HMM [-fh] -i <infile> -a <path | probability> -b <alignment>`

* `-f`: Add this flag to enable forward mode.
* `-h`: Print help information.
* `infile`: The input fastq file, in which every 2 sequences are considered as a pair.
* `path | probability`: The output path of inferred hidden states or probability in forward mode.
* `alignment`: The output path of alignment result, and is ignored in forward mode.

## Compilation
This project is managed by CMake, be sure you have a version over 3.23.
1. Go to the source code directory
2. `mkdir make-dir`
3. `cd make-dir`
4. `cmake ..`
5. `make`

Then find the executable file called "HMM" in make-dir.