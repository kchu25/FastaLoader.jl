mutable struct FASTA_DNA{S <: Real}
    N::Int
    L::Int
    acgt_freq::Vector{S}
    markov_bg_mat::Matrix{S}
    raw_data::Vector{String}
    data_matrix::Union{Array{S,3}, Array{S,2}}
    data_matrix_gpu::Union{CuArray{S,3}, CuArray{S,2}}
    data_matrix_bg::Union{Array{S,3}, Array{S,2}}
    labels::Union{Nothing, Vector{String}, Vector{Int}}

    function FASTA_DNA{S}(fasta_location::String; 
                        max_entries=max_num_read_fasta,
                        ryan_w_labels=false
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
    data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat = get_data_matrices(dna_read; FloatType=S);
    N = length(dna_read); L = Int(size(data_matrix,1)/4);
    data_matrix = reshape(data_matrix, 4*L, 1, N);
    new(        
        N,
        L,
        acgt_freq,
        markov_bg_mat,
        dna_read,
        data_matrix,
        cu(data_matrix),
        reshape(data_matrix_bg, 4*L, 1, N),
        labels
        )    
    end
end


get_N(d::FASTA_DNA) = d.N;
get_L(d::FASTA_DNA) = d.L;
get_data_matrix(d::FASTA_DNA) = d.data_matrix;
get_data_matrix_bg(d::FASTA_DNA) = d.data_matrix_bg;


