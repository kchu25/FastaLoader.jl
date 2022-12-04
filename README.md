# FastaLoader

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/FastaLoader.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/FastaLoader.jl/dev)
[![Build Status](https://github.com/kchu25/FastaLoader.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/FastaLoader.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/FastaLoader.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/FastaLoader.jl)

 
This is a package that provides subroutines that loads the DNA sequences in the specified fasta file. The DNA sequences are then transformed into some other useful information, e.g. one-hot/WYK encoded vectors, kmer-frequency preserved shuffled sequences, Markov background estimates, partitioned datasets for K-fold cross-validations (for fasta with labels), etc. for downstream machine learning tasks. As of now, we require all sequences in the fasta file to be the same length, and strings must be defined on DNA alphabets `{A,C,G,T}`.

# Usage

Coming Soon

<!-- 
## Loading fasta file for unsupervised learning tasks

Let's use JASPAR dataset [MA0463.1](https://jaspar.genereg.net/matrix/MA0463.1/) as an example; the fasta file that we can direcly download is [MA0463.1.sites](https://jaspar.genereg.net/download/data/2022/sites/MA0463.1.sites).

Once downloaded, you can load the this Fasta file by simply calling 
```
using FastaLoader

ma0463_1 = FASTA_DNA{Float32}(<path to MA0463.1.sites>)
``` -->




