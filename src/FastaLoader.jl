module FastaLoader

using SeqShuffle, CUDA

export FASTA_DNA
include("constants.jl")
include("helpers.jl")
include("fasta.jl")

end
