



# fasta = make_FASTA_DNA_w_splits(fp2);


# FastaLoader.get_train_fold_for_flux(fasta, 1) |> typeof
# q1,q2,q3,q4=FastaLoader.get_train_fold_for_flux(fasta, 1) 



# q2=get_test_set_ind(mcs);


# # this is for ryan's data now 
# mutable struct FASTA_DNA_folds{S <: Real}
#     N::Vector{Int}
#     L::Vector{Int}
#     raw_data::Vector{String}
#     data_matricies::Union{Vector{Array{S,3}}, Vector{Array{S,2}}}
#     data_matricies_gpu::Union{Vector{CuArray{S,3}}, Vector{CuArray{S,2}}}
#     data_matricies_bg::Union{Vector{Array{S,3}}, Vector{Array{S,2}}}
#     labels::Union{Nothing, Vector{Vector{String}}, Vector{Vector{Int}}}

#     function FASTA_DNA_folds{S}(fasta_location::String; fold=5)
#         dna_read = nothing; labels = nothing;
#         all_labels, all_dna_read = FastaLoader.reading(fasta_location; get_header=true, ryan_data=true);
#         dcount = FastaLoader.get_count_map(all_labels); # count the number data point asscociated with each label
#         ks = [k for k in keys(dcount)]; vals = [v for v in values(dcount)];
#         valid_labels = ks[vals .> label_count_thresh]; # labels that got more data to be considered valid
#         class_indicators = reduce(hcat, [all_labels .== v for v in valid_labels]); # bitarrays for class indicates
#         class_indices = [findall(@view class_indicators[:,i]) for i = 1:size(class_indicators,2)]
#         # shuffle 
#         shuffles_class_indices = shuffle.(class_indices);
#         splits = train_test_split.(shuffles_class_indices);
#         shuffles_class_indices_chunk = chunk.(shuffles_class_indices, fold)

#         # data_matrices = 
#     end
# end

# # println("hi")

# function split_train_test(header_labels, dna_reads)

# end
