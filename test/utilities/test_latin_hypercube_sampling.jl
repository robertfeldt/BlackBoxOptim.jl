@testset "Latin hypercube sampling" begin
    @test_throws DimensionMismatch BlackBoxOptim.Utils.latin_hypercube_sampling(Float64[], [1.0], 1)
    @test_throws DimensionMismatch BlackBoxOptim.Utils.latin_hypercube_sampling([1.0, 2.0], [1.0], 1)
    @test_throws ArgumentError BlackBoxOptim.Utils.latin_hypercube_sampling([2.0], [1.0], 1)

    samples0 = BlackBoxOptim.Utils.latin_hypercube_sampling([2.0], [2.0], 1)
    @test samples0 isa Matrix{Float64}
    @test size(samples0) == (1, 1)
    @test samples0[1,1] == 2.0

    samples1 = BlackBoxOptim.Utils.latin_hypercube_sampling([0.0], [1.0], 1)
    # since the output of latin_hypercube_sampling() would be used as population matrix
    # it should be of Matrix type, not just AbstractMatrix
    @test samples1 isa Matrix{Float64}
    @test size(samples1) == (1, 1)
    @test 0.0 <= samples1[1,1] <= 1.0

    samples2 = BlackBoxOptim.Utils.latin_hypercube_sampling([0.0], [1.0], 2)
    @test samples2 isa Matrix{Float64}
    @test size(samples2) == (1, 2)
    samples2 = sort(samples2, dims=2)
    @test 0.0 <= samples2[1,1] <= 0.5
    @test 0.5 <= samples2[1,2] <= 1.0

    samples3 = BlackBoxOptim.Utils.latin_hypercube_sampling([0.0, 2.0], [1.0, 3.0], 4)
    @test samples3 isa Matrix{Float64}
    @test size(samples3) == (2, 4)
    samples3 = sort(samples3, dims=2)
    @test 0.0 <= samples3[1,1] <= 0.25
    @test 0.25 <= samples3[1,2] <= 0.5
    @test 0.5 <= samples3[1,3] <= 0.75
    @test 0.75 <= samples3[1,4] <= 1.0
    @test 2.0 <= samples3[2,1] <= 2.25
    @test 2.25 <= samples3[2,2] <= 2.5
    @test 2.5 <= samples3[2,3] <= 2.75
    @test 2.75 <= samples3[2,4] <= 3.0
end
