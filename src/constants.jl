######## fasta load settings: ##################
const max_num_read_fasta = 100000;
const dat_t = Float32;          # data_matrix_type
const int_t = Int32;            # integer type

#= each label must associated with at least 200 sequences 
   for classification 
=#
const label_count_thresh = 250; 
