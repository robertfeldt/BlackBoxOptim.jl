NumTestRepetitions = 100

@testset "Adaptive differential evolution optimizer" begin

ss = RectSearchSpace(1, (0.0, 10.0))
fake_problem = FunctionBasedProblem(x -> 0.0, "test_problem", MinimizingFitnessScheme, ss)

ade = adaptive_de_rand_1_bin(fake_problem, ParamsDict(
    :Population => rand(1, 100)))

@testset "parameters adjust!()" begin
    for i in 1:NumTestRepetitions
        cur_cr, cur_f = BlackBoxOptim.crossover_parameters(ade.modify, ade.population, rand(1:popsize(ade.population)))
        @test (0.0 <= cur_cr <= 1.0)
        @test (0.0 <= cur_f <= 1.0)
        BlackBoxOptim.adjust!(ade.modify, 0, 1, 0.0, 0.0, false)
        # FIXME this fails too often; needs adjusting the distribution?
        #new_cr, new_f = BlackBoxOptim.crossover_parameters(ade.params, 1)
        #@fact new_cr != cur_cr --> true
        #@fact new_f != cur_f --> true
    end
end

# FIXME actually this test is not required as the standard DE already tests for that
@testset "ask()" begin
    for i in 1:NumTestRepetitions
        res = BlackBoxOptim.ask(ade)
        @test length(res) == 2

        trial, target = res

        @test ndims(trial.params) == 1
        @test (1 <= trial.index <= popsize(ade))
        @test in(trial.params, search_space(ade.embed))

        @test ndims(target.params) == 1
        @test (1 <= target.index <= popsize(ade))
        @test in(target.params, search_space(ade.embed))

        @test trial.index == target.index
    end
end

end
