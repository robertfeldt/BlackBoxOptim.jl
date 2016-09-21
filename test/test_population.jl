@testset "Population" begin

    @testset "FitPopulation" begin
        fs = MinimizingFitnessScheme
        p1 = FitPopulation(fs, 10, 2)
        @test isa(p1, FitPopulation)

        @test popsize(p1) == 10
        @test numdims(p1) == 2

        @test isnafitness(fitness(p1, 1), fs)
        @test isnafitness(fitness(p1, 4), fs)

        @testset "candidates pool" begin
            @test BlackBoxOptim.candi_pool_size(p1) == 0
            candi1 = BlackBoxOptim.acquire_candi(p1, 1)
            @test BlackBoxOptim.candi_pool_size(p1) == 0
            @test candi1.index == 1
            @test isnafitness(candi1.fitness, fs)

            candi2 = BlackBoxOptim.acquire_candi(p1, 2)
            @test BlackBoxOptim.candi_pool_size(p1) == 0
            @test candi2.index == 2
            @test isnafitness(candi2.fitness, fs)

            BlackBoxOptim.release_candi(p1, candi2)
            @test BlackBoxOptim.candi_pool_size(p1) == 1

            candi1.fitness = 5.0
            BlackBoxOptim.accept_candi!(p1, candi1)
            @test BlackBoxOptim.candi_pool_size(p1) == 2
            @test fitness(p1, 1) == 5.0
        end

        @testset "append!()" begin
            p2 = FitPopulation(fs, 5, 2)
            @test isa(p2, FitPopulation)

            candi22 = BlackBoxOptim.acquire_candi(p2, 2)
            candi22.fitness = 4.0
            BlackBoxOptim.accept_candi!(p2, candi22)

            append!(p1, p2)
            @test numdims(p1) == 2
            @test popsize(p1) == 15
            @test p1[11] == p2[1]
            @test fitness(p1, 12) == 4.0

            p3 = FitPopulation(fs, 5, 1)
            @test_throws DimensionMismatch append!(p1, p3)
        end
    end

end
