@testset "Latin hypercube sampling" begin
    samples = BlackBoxOptim.Utils.latin_hypercube_sampling([0.0], [1.0], 1)
    @test size(samples, 1) == 1
    @test (0.0 <= samples[1,1] <= 1.0)

    samples = BlackBoxOptim.Utils.latin_hypercube_sampling([0.0], [1.0], 2)
    @test size(samples, 2) == 2
    sorted = sort(samples, 2)
    @test (0.0 <= sorted[1,1] <= 0.5)
    @test (0.5 <= sorted[1,2] <= 1.0)

    samples = BlackBoxOptim.Utils.latin_hypercube_sampling([0.0, 2.0], [1.0, 3.0], 4)
    @test size(samples, 1) == 2
    @test size(samples, 2) == 4
    sorted = sort(samples[[1],:],2)
    @test (0.0 <= sorted[1,1] <= 0.25)
    @test (0.25 <= sorted[1,2] <= 0.5)
    @test (0.5 <= sorted[1,3] <= 0.75)
    @test (0.75 <= sorted[1,4] <= 1.0)
    s2 = sort(samples[[2],:],2)
    @test (2.0 <= s2[1,1] <= 2.25)
    @test (2.25 <= s2[1,2] <= 2.5)
    @test (2.5 <= s2[1,3] <= 2.75)
    @test (2.75 <= s2[1,4] <= 3.0)
end
