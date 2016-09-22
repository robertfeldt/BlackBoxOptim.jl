@testset "Crossover operators" begin

ss = symmetric_search_space(1, (0.0, 10.0))
fake_problem = FunctionBasedProblem(x -> 0.0, "test_problem", MinimizingFitnessScheme, ss)
DE = de_rand_1_bin(fake_problem, ParamsDict(
    :Population => collect(1.0:10.0)',
    :f => 0.4, :cr => 0.5, :NumParents => 3))

@testset "DiffEvoRandBin1" begin
    @testset "always copies from donor if length is 1" begin
        @test BlackBoxOptim.apply!(BlackBoxOptim.DiffEvoRandBin1(0.0, 0.0),
                                   [0.0], 4, DE.population, [1,2,3]) == [3.0]
    end

    @testset "always copies at least one element from donor" begin
        for i in 1:NumTestRepetitions
            len = rand(1:100)
            pop = rand(len,4)
            target = pop[:,1]
            saved_target = copy(target)
            res = BlackBoxOptim.apply!(BlackBoxOptim.DiffEvoRandBin1(0.0, 0.0),
                                       target, 1, pop, [2,3,4])
            @test (res === target)
            @test ndims( res ) == 1
            @test length( res ) == len
            @test sum(target .!= saved_target) == 1
        end
    end

    @testset "unlikely to copy everything if vectors are large" begin
        for i in 1:NumTestRepetitions
            len = 50 + rand(1:50)
            pop = rand(len,4)
            target = pop[:,1]
            saved_target = copy(target)
            res = BlackBoxOptim.apply!(BlackBoxOptim.DiffEvoRandBin1(0.1, 0.5),
                                       target, 1, pop, [2,3,4])
            @test any(target .!= saved_target)
            @test any(target .== saved_target)
        end
    end

    @testset "correctly modifies the parameters vector" begin
        @test BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                                                [0.0], 4, DE.population, [1,2,3]) == [3.0 + (0.4 * (1.0 - 2.0))]
        @test BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                                                [0.0], 4, DE.population, [2,3,1]) == [1.0 + (0.4 * (2.0 - 3.0))]
        @test BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                                                [0.0], 4, DE.population, [3,2,1]) == [1.0 + (0.4 * (3.0 - 2.0))]
        @test BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                                                [0.0], 4, DE.population, [1,3,5]) == [5.0 + (0.4 * (1.0 - 3.0))]
        @test BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                                                [0.0], 3, DE.population, [4,9,8]) == [8.0 + (0.4 * (4.0 - 9.0))]

        pop2 = reshape(collect(1.0:8.0), 4, 2)'
        @test BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.6),
                                                                [0.0,0.0], 4, pop2, [1,2,3]) == [3.0 + (0.6 * (1.0 - 2.0)),
                                                                                                 7.0 + (0.6 * (5.0 - 6.0))]
        @test BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.6),
                                                                [0.0,0.0], 4, pop2, [1,2,4]) == [4.0 + (0.6 * (1.0 - 2.0)),
                                                                                                 8.0 + (0.6 * (5.0 - 6.0))]
    end
end

@testset "MutationWrapper" begin
    ss = symmetric_search_space(2, (0.0, 10.0))
    pop = reshape(collect(1.0:8.0), 2, 4)
    gibbs = UniformMutation(ss)
    gibbs_wrapper = BlackBoxOptim.MutationWrapper(gibbs)
    @test numchildren(gibbs_wrapper) == 1
    @test numparents(gibbs_wrapper) == 1
    mut_res = BlackBoxOptim.apply!(gibbs_wrapper, [0.0, 0.0], 1, pop, [2])
    @test sum(mut_res .== BlackBoxOptim.viewer(pop, 2)) == 0
end

end
