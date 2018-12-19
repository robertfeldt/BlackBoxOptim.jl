@testset "Differential evolution optimizer" begin

ss = RectSearchSpace(1, (0.0, 10.0))
fake_problem = FunctionBasedProblem(x -> 0.0, "test_problem", MinimizingFitnessScheme, ss)
DE = de_rand_1_bin(fake_problem, ParamsDict(
    :Population => reshape(collect(1.0:10.0), (1, 10)),
    :f => 0.4, :cr => 0.5, :NumParents => 3))

@testset "ask()/tell!()" begin
    for i in 1:NumTestRepetitions
        res = BlackBoxOptim.ask(DE)
        @test BlackBoxOptim.candi_pool_size(BlackBoxOptim.population(DE)) == 0 # no candidates in the pool as we just exhausted it
        @test length(res) == 2

        trial, target = res

        @test ndims(trial.params) == 1
        @test (1 <= trial.index <= popsize(DE.population))
        @test in(trial.params, search_space(DE.embed))

        @test ndims(target.params) == 1
        @test (1 <= target.index <= popsize(DE.population))
        @test in(target.params, search_space(DE.embed))

        @test trial.index == target.index

        BlackBoxOptim.tell!(DE, res)
        @test BlackBoxOptim.candi_pool_size(BlackBoxOptim.population(DE)) == 2 # test that all candidates returned to the pool
    end
end

end
