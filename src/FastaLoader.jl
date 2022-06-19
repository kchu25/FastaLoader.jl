module FastaLoader

using SeqShuffle, CUDA, StatsBase, Random

export FASTA_DNA, 
       make_FASTA_DNA_w_splits,
       get_test_set_for_flux,
       get_train_fold_for_flux,
       fasta_reshape_for_flux!,
       fasta_reshape_to_orig!


include("constants.jl")
include("helpers.jl")
include("fasta.jl")
include("fasta_w_splits.jl")


end
