struct class_split
    train_folds::Vector{Vector{Int}}
    testset::Vector{Int}
end

struct multiple_class_splits
    splits::Vector{class_split}
    num_folds::Int # the number of folds in each train_folds in class_split in splits
end

mutable struct FASTA_DNA_w_splits{S <: Real}
    # mutable for matrix reshapings
    mcs::multiple_class_splits

    labels::Vector{String}
    label_indicators::BitMatrix # note that label indicators includes all labels and is not shuffled
    label_indicators_gpu::CuArray{S,2}
    raw_data::Vector{String}

    acgt_freq::Vector{S}
    markov_bg_mat::Matrix{S}

    data_matrix::Union{Array{S,3}, Array{S,2}}
    data_matrix_gpu::Union{CuArray{S,3}, CuArray{S,2}}
    data_matrix_bg::Union{Array{S,3}, Array{S,2}}
end

function equal_chunk(arr; num_chunks=n)
    # parition an array to n equal-sized chunks where the 
    # last chunk may be slightly larger
    len_arr = length(arr);
    chunk_size = div(len_arr, num_chunks)
    chunks = Vector{Vector{eltype(arr)}}();
    for i = 1:chunk_size:len_arr
        if i + 2*chunk_size - 1  > len_arr
            push!(chunks, arr[i:end]);
            break
        else
            push!(chunks, arr[i:i+chunk_size-1]);
        end
    end
    return chunks
end

function train_test_split(arr; split_ratio=0.85, folds=5)
    # return the data as k-fold cross-validation scheme
    # return (k-folds, test-set) as a class_split instance
    # a separate test-set is needed since we've used 
    # k-folds for hyperparameter selection
    train_set_size = Int(floor(split_ratio*length(arr)));
    return class_split(
                equal_chunk(arr[1:train_set_size]; num_chunks=folds), 
                arr[train_set_size+1:end]
                )
end

function get_ryan_fasta_str_labels(fasta_location::String)
    # return the labels and the corresponding dna strings in ryan's data
    all_labels, all_dna_read = FastaLoader.reading(fasta_location; get_header=true, ryan_data=true);
    return all_labels, all_dna_read
end

function read_ryan_fasta(all_labels::Vector{String}; 
                         shuffle_this=true, 
                         count_thresh=label_count_thresh)
    # Return a vector of shuffled array where each array
    # is the indices of the strings in fasta that correspond to each class,
    # and the one-hot encoded vector for the corresponding labels
    # This only returns classes that have a sufficient amount of data (count_thresh).
    dcount = FastaLoader.get_count_map(all_labels); # count the number data point asscociated with each label
    ks = [k for k in keys(dcount)]; vals = [v for v in values(dcount)];
    valid_labels = ks[vals .> count_thresh]; # labels that got more data to be considered valid    
    class_indicators = reduce(hcat, [all_labels .== v for v in valid_labels]); # bitarrays for class indicates
    class_indices = [findall(@view class_indicators[:,i]) for i = 1:size(class_indicators,2)]    
    shuffles_class_indices = shuffle_this ? shuffle.(class_indices) : class_indices;
    # note that class_indicators is not shuffled
    # TODO: remove allocations later
    return shuffles_class_indices, Array(class_indicators') 
end

function make_FASTA_DNA_w_splits(fp::String; 
                                 class_selector=read_ryan_fasta, 
                                 split_ratio=0.85,
                                 folds=5,
                                 flux=true,
                                 float_type=Float32)
    all_labels, all_dna_read = get_ryan_fasta_str_labels(fp);
    shuffles_class_indices, class_indicators = class_selector(all_labels);
    # split each class to have train and test set
    mcs = multiple_class_splits(
                train_test_split.(shuffles_class_indices; 
                                  split_ratio=split_ratio, 
                                  folds=folds),
                folds
                );
    data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat = FastaLoader.get_data_matrices(all_dna_read; 
                                                                                 FloatType=float_type);
    fws = FASTA_DNA_w_splits(mcs, 
                              all_labels,
                              class_indicators,
                              cu(class_indicators),
                              all_dna_read,
                              acgt_freq,
                              markov_bg_mat,
                              data_matrix,
                              cu(data_matrix),
                              data_matrix_bg
                              );
    flux && fasta_reshape_for_flux!(fws);
end

function get_test_set_ind(mcs::multiple_class_splits)
    test_set_ind = Vector{Int}();
    for class_split in mcs.splits
        append!(test_set_ind, class_split.testset)
    end
    return test_set_ind
end

function get_train_fold_ind(mcs::multiple_class_splits, fold::Int)
    @assert fold ≤ mcs.num_folds "The input integer fold must be smaller than pre-specified number of folds."
    valid_set_ind = Vector{Int}();
    train_set_ind = Vector{Int}();
    for class_split in mcs.splits # for each class
        for f = 1:mcs.num_folds
            if f != fold
                append!(train_set_ind, class_split.train_folds[f]);
            else
                append!(valid_set_ind, class_split.train_folds[f]);
            end
        end
    end
    return train_set_ind, valid_set_ind
end

function fasta_reshape_to_orig!(fws::FASTA_DNA_w_splits)
    L, _, N = size(fws.data_matrix);
    fws.data_matrix = reshape(fws.data_matrix, (L, N));
    fws.data_matrix_gpu = reshape(fws.data_matrix_gpu, (L, N));
    fws.data_matrix_bg = reshape(fws.data_matrix_bg, (L, N));
end

function fasta_reshape_for_flux!(fws::FASTA_DNA_w_splits)
    # reshape the data so that Flux dataloader can take it
    L, N = size(fws.data_matrix);
    fws.data_matrix = reshape(fws.data_matrix, (L,1,N));
    fws.data_matrix_gpu = reshape(fws.data_matrix_gpu, (L,1,N));
    fws.data_matrix_bg = reshape(fws.data_matrix_bg, (L,1,N));
end

function get_test_set_for_flux(fws::FASTA_DNA_w_splits; gpu=true)
    test_set_ind = get_test_set_ind(fws.mcs);
    if gpu
        return fws.data_matrix_gpu[:,:,test_set_ind], fws.label_indicators_gpu[:,test_set_ind]
    else
        return fws.data_matrix[:,:,test_set_ind], fws.label_indicators[:,test_set_ind]
    end
end

function get_train_fold_for_flux(fws::FASTA_DNA_w_splits, fold::Int; gpu=true)
    train_set_ind, valid_set_ind = get_train_fold_ind(fws.mcs, fold)
    if gpu
        return fws.data_matrix_gpu[:,:,train_set_ind],  fws.label_indicators_gpu[:,train_set_ind]
               fws.data_matrix_gpu[:,:,valid_set_ind],  fws.label_indicators_gpu[:,valid_set_ind]
    else
        return fws.data_matrix[:,:,train_set_ind], fws.label_indicators[:,train_set_ind]
               fws.data_matrix[:,:,valid_set_ind], fws.label_indicators[:,valid_set_ind]
    end
end

