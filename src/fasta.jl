const dna_meta_data = Vector{NamedTuple{(:str, :motif_where, :mode), 
                            Tuple{String, UnitRange{Int64}, Int64}}}


mutable struct FASTA_with_BigWig{S <: Real}
    N::Int                                              # number of dna strings in the training set
    L::Int                                              # length of the dna strings in the training set
    N_test::Int                                         # number of dna strings in the test set
    L_test::Int                                         # length of the dna strings in the test set
    raw_data_train::Vector{String}                      # raw data of the training set
    raw_data_test::Vector{String}                       # raw data of the test set

end


mutable struct FASTA_DNA_for_regressions{S <: Real}
    N::Int                                              # number of dna strings in the training set
    L::Int                                              # length of the dna strings in the training set
    N_test::Int                                         # number of dna strings in the test set
    L_test::Int                                         # length of the dna strings in the test set
    raw_data_train::Vector{String}                      # raw data of the training set
    raw_data_test::Vector{String}                       # raw data of the test set
    data_matrix::Union{Array{S,3}, Array{S,2}}          # data array (one-hot, cpu) of the training set
    data_matrix_test::Union{Array{S,3}, Array{S,2}}     # data array (one-hot, cpu) of the training set
    labels::Vector{S}                                   # training set string labels
    labels_test::Vector{S}                              # test set string labels

    function FASTA_DNA_for_regressions{S}(fasta::String; 
            train_test_ratio=0.8) where {S <: Real}

        # TODO: do cross-validation for the test set later
        n_train, n_test, labels_train, seqs_train, labels_test, seqs_test = 
            read_and_permute(fasta; 
                train_test_ratio=train_test_ratio, 
                parse_float_type=S)

        data_matrix = data_2_dummy(seqs_train; F=S);
        data_matrix_test = data_2_dummy(seqs_test; F=S);

        L = length(seqs_train[1]); 
        L_test = length(seqs_test[1]);

        new(
            n_train,
            L,
            n_test,
            L_test,
            seqs_train, 
            seqs_test, 
            reshape(data_matrix, (4*L, 1, n_train)),
            reshape(data_matrix_test, (4*L_test, 1, n_test)),
            labels_train,
            labels_test
        )
    end
end

mutable struct FASTA_DNA_for_classifications{S <: Real}
    N::Int                                              # number of dna strings in the training set
    L::Int                                              # length of the dna strings in the training set
    N_test::Int                                         # number of dna strings in the test set
    L_test::Int                                         # length of the dna strings in the test set
    num_classes::Int                                    # number of classes
    raw_data_train::Vector{String}                      # raw data of the training set
    raw_data_test::Vector{String}                       # raw data of the test set
    data_matrix::Union{Array{S,3}, Array{S,2}}          # data array (one-hot, cpu) of the training set
    data_matrix_test::Union{Array{S,3}, Array{S,2}}     # data array (one-hot, cpu) of the training set
    labels::BitMatrix                                   # training set string labels
    labels_test::BitMatrix                              # test set string labels

    # TODO: labels to one hot encoding

    function FASTA_DNA_for_classifications{S}(fasta_train::String, fasta_test::String) where {S <: Real}
        labels, raw_data_train = reading_for_DNA_classification(fasta_train)
        labels_test, raw_data_test = reading_for_DNA_classification(fasta_test)    
        data_matrix = data_2_dummy(raw_data_train; F=S);
        data_matrix_test = data_2_dummy(raw_data_test; F=S);
        N = size(labels, 2); L = length(raw_data_train[1]); 
        N_test = size(labels_test, 2); L_test = length(raw_data_test[1]);
        new(
            N,
            L,
            N_test,
            L_test,
            length(unique(labels)),
            raw_data_train, 
            raw_data_test, 
            reshape(data_matrix, (4*L, 1, N)),
            reshape(data_matrix_test, (4*L_test, 1, N_test)),
            labels,
            labels_test
        )
    end
end

mutable struct FASTA_DNA{S <: Real}
    N::Int
    L::Int
    acgt_freq::Vector{S}
    markov_bg_mat::Matrix{S}
    raw_data::Vector{String}
    raw_data_test::Vector{String}
    data_matrix::Union{Array{S,3}, Array{S,2}}
    data_matrix_gpu::Union{CuArray{S,3}, CuArray{S,2}, Nothing}
    data_matrix_bg::Union{Array{S,3}, Array{S,2}}
    data_matrix_bg_gpu::Union{CuArray{S,3}, CuArray{S,2}, Nothing}
    labels::Union{Nothing, Vector{String}, Vector{Int}}
    meta_data::Union{Nothing, dna_meta_data}
    acgt_freq_test::Union{Nothing, Vector{S}}
    markov_bg_mat_test::Union{Nothing, Matrix{S}}
    data_matrix_test::Union{Nothing, Array{S,3}, Array{S,2}}
    data_matrix_bg_test::Union{Nothing, Array{S,3}, Array{S,2}}
    N_test::Int

    function FASTA_DNA{S}(fasta_location::String; 
                        max_entries=max_num_read_fasta,
                        k_train=1, k_test=2, # kmer frequency in the test set 
                        train_test_split_ratio=0.9,
                        shuffle=true
                        ) where {S <: Real}       
        dna_read = nothing; labels = nothing;
        dna_read = read_fasta(fasta_location; max_entries);
        data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat,
            data_matrix_test, data_matrix_bg_test, _, acgt_freq_test, 
                markov_bg_mat_test, N_train, N_test, train_set_inds, test_set_inds = 
                get_data_matrices(dna_read; k_train=k_train, k_test=k_test, 
                                  train_test_split_ratio=train_test_split_ratio, 
                                  shuffle=shuffle, 
                                  FloatType=S);
        L = Int(size(data_matrix,1)/4);
        data_matrix = reshape(data_matrix, 4*L, 1, N_train);
        data_matrix_test = reshape(data_matrix_test, 4*L, 1, N_test)
        data_matrix_bg = reshape(data_matrix_bg, 4*L, 1, N_train)
        new(        
            N_train,
            L,
            acgt_freq,
            markov_bg_mat,
            dna_read[train_set_inds],
            dna_read[test_set_inds],
            data_matrix,
            nothing,
            data_matrix_bg,
            nothing,
            labels,
            nothing,
            acgt_freq_test,
            markov_bg_mat_test,
            data_matrix_test,
            data_matrix_bg_test,
            N_test
            )
    end
end

mutable struct FASTA_DNA_JASPAR{S <: Real}
    N::Int
    L::Int
    acgt_freq::Vector{S}
    markov_bg_mat::Matrix{S}
    raw_strs::Vector{String}
    raw_strs_test::Vector{String}
    data_matrix::Union{Array{S,3}, Array{S,2}}
    data_matrix_gpu::Union{CuArray{S,3}, CuArray{S,2}}
    data_matrix_bg::Union{Array{S,3}, Array{S,2}}
    data_matrix_bg_gpu::Union{CuArray{S,3}, CuArray{S,2}}
    raw_data::Union{Nothing, dna_meta_data}
    raw_data_test::Union{Nothing, dna_meta_data}
    acgt_freq_test::Union{Nothing, Vector{S}}
    markov_bg_mat_test::Union{Nothing, Matrix{S}}
    data_matrix_test::Union{Nothing, Array{S,3}, Array{S,2}}
    data_matrix_bg_test::Union{Nothing, Array{S,3}, Array{S,2}}
    N_test::Int

    # constructor for JASPAR datasets
    function FASTA_DNA_JASPAR{S}(filepath::String; 
                                 max_entries=max_num_read_fasta,
                                 k_train=1, k_test=2, # kmer frequency in the test set
                                 train_test_split_ratio=0.9,
                                 shuffle=true
                                 ) where {S <: Real}
        dna_reads = reading(filepath; max_entries);   
        raw_data = dna_meta_data()
        num_ground_truth = 0;

        for str in dna_reads
            characters = [i for i in str];
            motif_found = findall(isuppercase.(characters) .== 1);
            num_ground_truth += isempty(motif_found) ? 0 : 1;
            motif_start = minimum(motif_found); motif_end = maximum(motif_found);
            motif_range =  motif_start:motif_end;        
            push!(raw_data, (str=str, motif_where=motif_range, mode=1));
        end

        dna_reads = uppercase.(dna_reads);

        data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat,
            data_matrix_test, data_matrix_bg_test, _, acgt_freq_test, 
                markov_bg_mat_test, N_train, N_test, train_set_inds, test_set_inds = 
                get_data_matrices(dna_reads; k_train=k_train, k_test=k_test,
                                train_test_split_ratio=train_test_split_ratio, 
                                FloatType=S,
                                shuffle=shuffle
                                );
        L = Int(size(data_matrix,1)/4);
        data_matrix = reshape(data_matrix, 4*L, 1, N_train);
        data_matrix_test = reshape(data_matrix_test, 4*L, 1, N_test)
        data_matrix_bg = reshape(data_matrix_bg, 4*L, N_train)
        new(N_train, 
            L, 
            acgt_freq, 
            markov_bg_mat, 
            dna_reads[train_set_inds], 
            dna_reads[test_set_inds],
            data_matrix, 
            cu(data_matrix), 
            data_matrix_bg,
            cu(data_matrix_bg),
            raw_data[train_set_inds],
            raw_data[test_set_inds],
            acgt_freq_test,
            markov_bg_mat_test,
            data_matrix_test,
            data_matrix_bg_test,
            N_test
           )
    end
end

get_N(d::FASTA_DNA) = d.N_train;
get_L(d::FASTA_DNA) = d.L;
get_data_matrix(d::FASTA_DNA) = d.data_matrix;
get_data_matrix_bg(d::FASTA_DNA) = d.data_matrix_bg;




