struct FASTA_DNA{S <: Real}
    N::Int
    L::Int
    acgt_freq::Vector{S}
    markov_bg_mat::Matrix{S}
    raw_data::Vector{String}
    data_matrix::Array{S,4}
    data_matrix_bg::Array{S,4}
    # data_matrix_bg_prob::Union{Nothing, Matrix{S}}

    function FASTA_DNA{S}(fasta_location::String, 
                        max_entries=max_num_read_fasta
                        ) where {S <: Real}       
    dna_read = read_fasta(fasta_location; max_entries);
    data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat = get_data_matrices(dna_read; FloatType=S);
    N = length(dna_read); L = Int(size(data_matrix,1)/4);
    new(        
        N,
        L,
        acgt_freq,
        markov_bg_mat,
        dna_read,
        reshape(data_matrix, 4*L, 1, 1, N),
        reshape(data_matrix_bg, 4*L, 1, 1, N))    
    end
end

get_N(d::FASTA_DNA) = d.N;
get_L(d::FASTA_DNA) = d.L;
get_data_matrix(d::FASTA_DNA) = d.data_matrix;
get_data_matrix_bg(d::FASTA_DNA) = d.data_matrix_bg;


