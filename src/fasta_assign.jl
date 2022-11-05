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
            this_read = join(uppercase.(splits[2:end]));  
            push!(dna_reads, this_read);
        end
    end
    return header_tuples, dna_reads
end

function read_supp3(supp3::String, 
                    supp3_fasta::String; 
                    twoE_oneS=false, 
                    strong_vs_silence=true,
                    reverse_comp=true)
    # get the sequences from Ryan's dataset                    
    @assert sum([twoE_oneS, strong_vs_silence]) ≤ 1 "only one extraction condition can be true"

    # read the labels into dataframe
    df3 = DataFrame(CSV.File(supp3)); 
    classes_we_want = nothing;
    group_name_WT_binary = nothing;

    if twoE_oneS
        group_name_WT_binary = map(x-> x=="Weak enhancer" || x == "Strong enhancer" ? "Enhancer" : x, df3.group_name_WT)
        classes_we_want = group_name_WT_binary .∈ [["Enhancer", "Silencer"]]
    elseif strong_vs_silence
        classes_we_want = df3.group_name_WT .∈ [["Strong enhancer", "Silencer"]]
    else
        classes_we_want = df3.group_name_WT .∈ [["Weak enhancer", "Strong enhancer",  "Silencer"]];
    end
   
    # get the DNA sequences ####################
    _, d = reading_w_header(supp3_fasta);
    ############################################

    # setup the labels and sequences
    labels = nothing; seqs = nothing;
    if twoE_oneS
        labels, seqs = String.(group_name_WT_binary[classes_we_want]), d[classes_we_want];
    else
        labels, seqs = String.(df3.group_name_WT[classes_we_want]), d[classes_we_want];
    end

    if reverse_comp
        labels = [labels..., labels...] 
        seqs = [seqs..., reverse_complement.(seqs)...];
    end

    return labels, seqs
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
            this_read = join(uppercase.(splits[2:end]));        
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