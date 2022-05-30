struct FASTA_DNA{T <: Integer, S <: Real}
    N::T
    L::T
    acgt_freq::Vector{S}
    markov_bg_mat::Matrix{S}
    raw_data::Vector{String}
    data_matrix::Matrix{S}
    data_matrix_bg::Matrix{S}
    data_matrix_bg_prob::Union{Nothing, Matrix{S}}
    target_folder::String

    function FASTA_DNA{T, S}(fasta_location::String, 
                        max_entries=max_num_read_fasta
                        ) where {T <: Integer, S <: Real}       
    dna_read = read_fasta(fasta_location; max_entries);
    data_matrix, data_matrix_bg, data_bg_prob, acgt_freq, markov_bg_mat = get_data_matrices(dna_read; FloatType=S);
    new(        
        T(length(dna_read)),
        T(size(data_matrix,1)/4),
        acgt_freq,
        markov_bg_mat,
        dna_read,
        data_matrix,
        data_matrix_bg,
        data_bg_prob)    
    end
end

get_N(d::FASTA_DNA) = d.N;
get_L(d::FASTA_DNA) = d.L;
get_data_matrix(d::FASTA_DNA) = d.data_matrix;
get_data_matrix_bg(d::FASTA_DNA) = d.data_matrix_bg;
findparam_real(d::FASTA_DNA{T,S}) where {T,S} = S;


