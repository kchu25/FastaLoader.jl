const dna_meta_data = Vector{NamedTuple{(:str, :motif_where, :mode), 
                            Tuple{String, UnitRange{Int64}, Int64}}}

mutable struct FASTA_DNA{S <: Real}
    N::Int
    L::Int
    acgt_freq::Vector{S}
    markov_bg_mat::Matrix{S}
    raw_data::Vector{String}
    raw_data_test::Vector{String}
    data_matrix::Union{Array{S,3}, Array{S,2}}
    data_matrix_gpu::Union{CuArray{S,3}, CuArray{S,2}}
    data_matrix_bg::Union{Array{S,3}, Array{S,2}}
    data_matrix_bg_gpu::Union{CuArray{S,3}, CuArray{S,2}}
    labels::Union{Nothing, Vector{String}, Vector{Int}}
    meta_data::Union{Nothing, dna_meta_data}
    acgt_freq_test::Union{Nothing, Vector{S}}
    markov_bg_mat_test::Union{Nothing, Matrix{S}}
    data_matrix_test::Union{Nothing, Array{S,3}, Array{S,2}}
    data_matrix_bg_test::Union{Nothing, Array{S,3}, Array{S,2}}
    N_test::Int

    function FASTA_DNA{S}(fasta_location::String; 
                        max_entries=max_num_read_fasta,
                        ryan_w_labels=false,
                        k_train=1, k_test=2, # kmer frequency in the test set 
                        train_test_split_ratio=0.9,
                        shuffle=true
                        ) where {S <: Real}       
        dna_read = nothing; labels = nothing;
        if ryan_w_labels
            all_labels, all_dna_read = reading(fasta_location; get_header=true, ryan_data=true);
            dcount = get_count_map(all_labels); # count the number data point asscociated with each label
            ks = [k for k in keys(dcount)]; vals = [v for v in values(dcount)];
            valid_labels = ks[vals .> label_count_thresh]; # labels that got more data to be considered valid
            indicators = map(x-> x ∈ valid_labels ? true : false, all_labels);
            labels = all_labels[indicators];
            dna_read = [uppercase(i) for i in all_dna_read[indicators]];
        else
            dna_read = read_fasta(fasta_location; max_entries);
        end
        data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat,
            data_matrix_test, data_matrix_bg_test, _, acgt_freq_test, 
                markov_bg_mat_test, N_train, N_test, train_set_inds, test_set_inds = 
                get_data_matrices(dna_read; k_train=k_train, k_test=k_test, 
                                  train_test_split_ratio=train_test_split_ratio, 
                                  shuffle=shuffle, 
                                  FloatType=S);
        # data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat = get_data_matrices(dna_read; FloatType=S);
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
            cu(data_matrix),
            data_matrix_bg,
            cu(data_matrix_bg),
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




