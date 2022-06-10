function reading(filepath::String, max_entries=max_num_read_fasta)
    # read the file
    reads = read(filepath, String);
    # process the fasta file to get the DNA part
    # rid of sequences that contains "n"
    dna_reads = Vector{String}();
    for i in split(reads, '>')
        if !isempty(i)
            this_read = join(split(i, "\n")[2:end]);        
            if !occursin("N", this_read) && !occursin("n", this_read)
                push!(dna_reads, this_read);
            end
        end
    end    
    # dna_reads = [join(split(i, "\n")[2:end]) for i in split(reads, '>') if !isempty(i) && !occursin("N", i)]; 
    dna_reads = length(dna_reads) > max_entries ? dna_reads[1:max_entries] : dna_reads;
    return dna_reads
end

function read_fasta(filepath::String; 
                    max_entries=max_num_read_fasta
                    )::Vector{String}
    #= read a fasta file =#
    dna_reads = reading(filepath, max_entries);   
    # convert all DNA seqeunce to uppercase
    return [uppercase(i) for i in dna_reads]
end

function dna2dummy(dna_string::String, dummy::Dict; F=Float32)    
    v = Array{F,1}(undef, 4*length(dna_string));
    for (index, alphabet) in enumerate(dna_string)
        start = (index-1)*4+1;
        v[start:start+3] = dummy[uppercase(alphabet)];
    end
    return v
end

#=
get the set of dummy-vectors from a set of dna-strings
the dummy-vectors are all of same length (for now)
=#
function data_2_dummy(dna_strings::Vector{String}; F=Float32)

    dummy = Dict('A'=>Array{F}([1, 0, 0, 0]), 
                   'C'=>Array{F}([0, 1, 0, 0]),
                   'G'=>Array{F}([0, 0, 1, 0]), 
                   'T'=>Array{F}([0, 0, 0, 1]));

    how_many_strings = length(dna_strings);
    @assert how_many_strings != 0 "There aren't DNA strings found in the input";
    _len_ = length(dna_strings[1]); # length of each dna string in data    
    _S_ = Array{F, 2}(undef, (4*_len_, how_many_strings));
    for i = 1:how_many_strings
        length(dna_strings[i]) == _len_ && (@inbounds _S_[:, i] = dna2dummy(dna_strings[i], dummy))
    end
    return _S_
end

function get_data_matrices(dna_read; k=2, FloatType=dat_t)
    shuffled_dna_read = seq_shuffle.(dna_read; k=k);
    data_matrix = data_2_dummy(dna_read; F=FloatType);
    data_matrix_bg = data_2_dummy(shuffled_dna_read; F=FloatType);

    # estimate the Markov background (order 1)
    acgt_freq, markov_mat = est_1st_order_markov_bg(shuffled_dna_read; F=FloatType);
    data_bg_prob = SeqShuffle.assign_bg_prob(shuffled_dna_read, markov_mat, acgt_freq);

    return data_matrix, data_matrix_bg, data_bg_prob, acgt_freq, markov_mat
end

# println("!23")