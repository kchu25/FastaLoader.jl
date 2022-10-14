module FastaLoader

using SeqShuffle, CUDA, StatsBase, Random, DataFrames, CSV

export FASTA_DNA, 
       FASTA_DNA_for_classifications, 
       make_FASTA_DNA_w_splits,
       make_FASTA_DNA_w_splits_activity,
       get_test_set_for_flux,
       get_train_fold_for_flux,
       fasta_reshape_for_flux!,
       fasta_reshape_to_orig!,
       reading,
       FASTA_DNA_JASPAR


include("constants.jl")
include("helpers.jl")
include("fasta.jl")
include("fasta_assign.jl")
include("fasta_w_splits.jl")


end
