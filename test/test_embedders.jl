@testset "Embedding operators" begin

@testset "RandomBound" begin
    @testset "does nothing if within bounds" begin
        @test BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                                                            [0.0], [0.0]) == [0.0]

        @test BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0), (10.0, 15.0)]),
                                                            [0.0, 11.4], [0.1, 12.3] ) == [0.0, 11.4]
    end

    @testset "bounds if lower than min bound" begin
        @test BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                                                            [-0.1], [0.0]) == [0.0]

        res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                                                            [-0.1], [0.5])
        @test res[1] >= 0.0
        @test res[1] <= 0.5

        res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(-11.0, 1.0), (0.0, 1.0)]),
                                                            [-11.1, 0.5], [-10.8, 0.5])
        @test (-10.8 >= res[1] >= -11.0)
        @test res[2] == 0.5

        res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(30.0, 60.0), (-102.0, -1.0)]),
                                                            [50.4, -103.1], [49.6, -101.4])
        @test res[1] == 50.4
        @test (-101.4 >= res[2] >= -102.0)
    end

    @testset "bounds if higher than max bound" begin
        @test BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                                                            [1.1], [1.0]) == [1.0]

        res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(-10.0, 96.0)]),
                                                            [97.0], [95.0])
        @test (95.0 <= res[1] <= 96.0)
    end
end

end
