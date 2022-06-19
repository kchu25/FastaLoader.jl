module FastaLoader

using SeqShuffle, CUDA, StatsBase, Random

export FASTA_DNA, make_FASTA_DNA_w_splits

include("constants.jl")
include("helpers.jl")
include("fasta.jl")
include("fasta_w_splits.jl")


end
