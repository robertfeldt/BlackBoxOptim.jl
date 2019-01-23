@testset "Population" begin

    @testset "FitPopulation" begin
        fs = MinimizingFitnessScheme
        p1 = FitPopulation(fs, 10, 2)
        @test isa(p1, FitPopulation)

        @test popsize(p1) == 10
        @test numdims(p1) == 2

        @test isnafitness(fitness(p1, 1), fs)
        @test isnafitness(fitness(p1, 4), fs)

        @testset "accessing individuals" begin
            @test isa(p1[1], Vector{Float64})
            @test length(p1[1]) == numdims(p1)
            @test isa(BlackBoxOptim.viewer(p1, 1), AbstractVector{Float64})
            @test length(BlackBoxOptim.viewer(p1, 1)) == numdims(p1)

            @test isa(p1[popsize(p1)], Vector{Float64}) # last solution vector
            @test_throws BoundsError p1[0]
            @test_throws BoundsError p1[popsize(p1)+1]
            @test_throws BoundsError BlackBoxOptim.viewer(p1, 0)
            @test_throws BoundsError BlackBoxOptim.viewer(p1, popsize(p1)+1)
            rand_solution_idx = rand(2:(popsize(p1)-1))
            @test isa(p1[rand_solution_idx], Array{Float64, 1}) # random solution vector
        end

        @testset "accessing individuals fitness" begin
            # and to access their fitness values:
            @test isa(fitness(p1, 1), Float64)
            @test isa(fitness(p1, popsize(p1)), Float64)
            rand_solution_idx = rand(2:(popsize(p1)-1))
            @test_throws BoundsError fitness(p1, 0)
            @test_throws BoundsError fitness(p1, popsize(p1)+1)
        end

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
