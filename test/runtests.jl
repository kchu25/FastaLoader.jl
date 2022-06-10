using FastaLoader
using Test

@testset "FastaLoader.jl" begin
    fp = "MA0463.1.sites";
    @test isfile(fp)

    # data = FASTA_DNA{Float32}(fp); # github doesn't have CUDA
    # println("Data $fp loaded."); 
    
    # @test data.L == 114
end
