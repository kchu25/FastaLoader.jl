# TODO add read with complement
function reading_w_header(filepath::String)
    reads = read(filepath, String);
    dna_reads = Vector{String}();
    header_tuples = Vector{Tuple{String,String,String}}();
    for i in split(reads, '>')
        if !isempty(i)
            splits = split(i, "\n");
            header_split = split(splits[1],"-");
            e = String(header_split[2]);
            chr, s = split(header_split[1],":");
            chr, s = String(chr), String(s);
            push!(header_tuples, (chr,s,e));
            this_read = join(splits[2:end]);  
            push!(dna_reads, this_read);
        end
    end
    return header_tuples, dna_reads
end


function reading_w_chr_loc(filepath::String;
                 max_entries=max_num_read_fasta,
                 get_header=false,
                 ryan_data=false)
    # read the file
    reads = read(filepath, String);
    # process the fasta file to get the DNA part
    # rid of sequences that contains "n"
    dna_reads = Vector{String}();
    header_tuples = Vector{Tuple{String,String,String}}();
    header_labels = Vector{String}();
    for i in split(reads, '>')
        if !isempty(i)
            splits = split(i, "\n");
            header_split = split(splits[1],"_");
            chr, s, e = split(header_split[1],"-")
            chr, s, e = String(chr), String(s), String(e);
            this_read = join(splits[2:end]);        
            if !occursin("N", this_read) && !occursin("n", this_read)
                push!(dna_reads, this_read);
                ryan_data && push!(header_labels, header_split[2])
                ryan_data && push!(header_tuples, (chr, s, e))
            end
        end
    end    
    # dna_reads = [join(split(i, "\n")[2:end]) for i in split(reads, '>') if !isempty(i) && !occursin("N", i)]; 
    dna_reads = length(dna_reads) > max_entries ? dna_reads[1:max_entries] : dna_reads;
    
    if get_header && ryan_data
        @assert length(dna_reads) == length(header_labels)
        return header_labels, header_tuples,  dna_reads
    else
        return dna_reads
    end
end