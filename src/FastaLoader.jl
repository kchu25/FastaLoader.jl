module FastaLoader

using SeqShuffle, CUDA, StatsBase, Random, DataFrames, CSV

export FASTA_DNA, 
       FASTA_DNA_for_regressions,
       FASTA_DNA_regression2,
       FASTA_DNA_for_classifications, 
       get_test_set_for_flux,
       get_train_fold_for_flux,
       fasta_reshape_for_flux!,
       fasta_reshape_to_orig!,
       reading,
       FASTA_DNA_JASPAR

include("constants.jl")
include("helpers.jl")
include("fasta.jl")


end
