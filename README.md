# HMM Viterbi Sequence Aligner

![cmake workflow](https://github.com/wykswr/HMM/actions/workflows/cmake.yml/badge.svg)

## Algorithm

This model utilizes paired hidden markov model (HMM) to find the alignment between 2 sequences.

### Parameters

#### Transition matrix

|     | 0   | 1              | 2               | 3               | 4              |
|-----|-----|----------------|-----------------|-----------------|----------------|
| 0   | $0$ | $\frac{1}{3}$  | $\frac{1}{3}$   | $\frac{1}{3}$   | $0$            |
| 1   | $0$ | $\frac{7}{10}$ | $\frac{1}{10}$  | $\frac{1}{10}$  | $\frac{1}{10}$ |
| 2   | $0$ | $\frac{1}{4}$  | $\frac{13}{20}$ | $0$             | $\frac{1}{10}$ |
| 3   | $0$ | $\frac{1}{4}$  | $0$             | $\frac{13}{20}$ | $\frac{1}{10}$ |
| 4   | $0$ | $0$            | $0$             | $0$             | $0$            |

#### $P_{seq_i, seq_j}$ and $P_{seq_i}$

$seq_i \quad ==\quad seq_j,\quad P_{seq_i, seq_j} = \frac{4}{28}$

$seq_i \quad != \quad seq_j,\quad P_{seq_i, seq_j} = \frac{1}{28}$

$\forall seq_i,\quad P_{seq_i} = \frac{1}{4}$

## Usage

The model is implemented as a UNIX command line tool, with syntax:

`HMM [-fh] -i <infile> -a <path | probability> -b <alignment>`

* `-f`: Add this flag to enable forward mode.
* `-h`: Print help information.
* `infile`: The input fastq file, in which every 2 sequences are considered as a pair.
* `path | probability`: The output path of inferred hidden states or probability in forward mode.
* `alignment`: The output path of alignment result, and is ignored in forward mode.

For instance, put pair.txt in the same directory of executable file called "HMM":
* Get hidden chain and alignment: `./HMM -i ./pair.txt -a ./path_alignment_model.txt -b ./predicted_alignment.txt`
* Get forward probability: `./HMM -i ./pair.txt -a ./prob_alignment_model.txt -f`
## Compilation

This project is managed by CMake, be sure you have a version over 3.23.

1. Go to the source code directory
2. `mkdir make-dir`
3. `cd make-dir`
4. `cmake ..`
5. `make`

Then find the executable file called "HMM" in make-dir.