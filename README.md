# FastaLoader

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/FastaLoader.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/FastaLoader.jl/dev)
[![Build Status](https://github.com/kchu25/FastaLoader.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/FastaLoader.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/FastaLoader.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/FastaLoader.jl)

 
This is a package that provides subroutines that loads the DNA sequences in the specified fasta file. The DNA sequences are then transformed into some other useful information, e.g. one-hot/WYK encoded vectors, kmer-frequency preserved shuffled sequences, Markov background estimates, partitioned datasets for K-fold cross-validations, etc. for downstream machine learning tasks. As of now, we require all sequences in the fasta file to be the same length, and strings must be defined on DNA alphabets `{A,C,G,T}`.

# Usage

Coming soon