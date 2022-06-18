using FastaLoader
using Test

@testset "FastaLoader.jl" begin
    fp = "MA0463.1.sites";
    fp2 = "ryan.txt";
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
end
