using FastaLoader
using Test

@testset "FastaLoader.jl" begin
    fp =  "MA0463.1.sites";
    @test isfile(fp)
    dna_reads = FastaLoader.read_fasta(fp);
    raw_data = FastaLoader.dna_meta_data()
    num_ground_truth = 0;

    data_matrix, data_matrix_bg, _, acgt_freq, markov_bg_mat,
    data_matrix_test, data_matrix_bg_test, _, acgt_freq_test, 
    markov_bg_mat_test, N_train, N_test, train_set_inds, test_set_inds = 
        FastaLoader.get_data_matrices(dna_reads; k_train=1, k_test=2, 
            train_test_split_ratio=0.9, shuffle=true, FloatType=Float32);

    @test typeof(dna_reads) == Vector{String}
    @test isapprox(sum(acgt_freq), 1; atol=0.001)
    @test isapprox(sum(acgt_freq_test), 1; atol=0.001)
    @test all(sum(markov_bg_mat, dims=2) .≈ 1)
    @test all(sum(markov_bg_mat_test, dims=2) .≈ 1)
end

