using FastaLoader
using Test

@testset "FastaLoader.jl" begin
    fp = "test/MA0463.1.sites";
    fp2 = "ryan.txt";
    fp3 = "supp3.txt"
    fp4 = "coor.fa"
    @test isfile(fp)
    @test isfile(fp2)

    # data = FASTA_DNA{Float32}(fp); # github doesn't have CUDA

    all_labels, all_dna_read = FastaLoader.reading(fp2; get_header=true, ryan_data=true);
    dcount = FastaLoader.get_count_map(all_labels); # count the number data point asscociated with each label
    ks = [k for k in keys(dcount)]; vals = [v for v in values(dcount)];
    valid_labels = ks[vals .> 200]; # labels that got more data to be considered valid
    indicators = map(x-> x ∈ valid_labels ? true : false, all_labels);
    labels = all_labels[indicators];
    dna_read = [uppercase(i) for i in all_dna_read[indicators]];

    @test typeof(dna_read) == Vector{String}
    @test length(labels) == length(dna_read)
    @test all(sort(unique(all_labels[indicators])) .== sort(valid_labels))

    split_ratio = 0.85; folds = 5;
    all_labels, all_dna_read = FastaLoader.get_ryan_fasta_str_labels(fp2);
    shuffles_class_indices, valid_labels, class_indicators = FastaLoader.read_ryan_fasta(all_labels);
    # split each class to have train and test set
    mcs = FastaLoader.multiple_class_splits(
                FastaLoader.train_test_split.(shuffles_class_indices; 
                                  split_ratio=split_ratio, 
                                  folds=folds),
                folds
                );
    test_set_ind = FastaLoader.get_test_set_ind(mcs);
    # these test that the indicator should always only contain valid labels 
    @test sum(sum(class_indicators[:,test_set_ind], dims=1) .== 0) == 0    
    for fold = 1:folds
        train_set_ind, valid_set_ind = FastaLoader.get_train_fold_ind(mcs, fold);
        @test sum(sum(class_indicators[:,train_set_ind], dims=1) .== 0) == 0
        @test sum(sum(class_indicators[:,valid_set_ind], dims=1) .== 0) == 0
    end

    # all_labels, all_dna_read = FastaLoader.read_supp3(fp3, fp4);
    # shuffles_class_indices, valid_labels, class_indicators = FastaLoader.read_ryan_fasta(all_labels);
    # mcs = FastaLoader.multiple_class_splits(
    #             FastaLoader.train_test_split.(shuffles_class_indices; 
    #                               split_ratio=split_ratio, 
    #                               folds=folds),
    #             folds
    #             );
    # test_set_ind = FastaLoader.get_test_set_ind(mcs);
    # @test sum(sum(class_indicators[:,test_set_ind], dims=1) .== 0) == 0    
    # for fold = 1:folds
    #     train_set_ind, valid_set_ind = FastaLoader.get_train_fold_ind(mcs, fold);
    #     @test sum(sum(class_indicators[:,train_set_ind], dims=1) .== 0) == 0
    #     @test sum(sum(class_indicators[:,valid_set_ind], dims=1) .== 0) == 0
    # end            
    
    # all_labels, all_dna_read = FastaLoader.read_supp3(fp3, fp4;  twoE_oneS=true, 
    #                                                     strong_vs_silence=false,
    #                                                     reverse_comp=true);
    # shuffles_class_indices, valid_labels, class_indicators = FastaLoader.read_ryan_fasta(all_labels);
    # mcs = FastaLoader.multiple_class_splits(
    #             FastaLoader.train_test_split.(shuffles_class_indices; 
    #                               split_ratio=split_ratio, 
    #                               folds=folds),
    #             folds
    #             );
    # test_set_ind = FastaLoader.get_test_set_ind(mcs);
    # @test sum(sum(class_indicators[:,test_set_ind], dims=1) .== 0) == 0    
    # for fold = 1:folds
    #     train_set_ind, valid_set_ind = FastaLoader.get_train_fold_ind(mcs, fold);
    #     @test sum(sum(class_indicators[:,train_set_ind], dims=1) .== 0) == 0
    #     @test sum(sum(class_indicators[:,valid_set_ind], dims=1) .== 0) == 0
    # end        

    # fp5 = "test/MA0599.1.sites"
    # dna_read = FastaLoader.read_fasta(fp5; max_entries=10000000);
    # data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat = FastaLoader.get_data_matrices(dna_read; FloatType=Float32);
    # my_read = read_fasta_that_has_ground_truth(fp);

end

