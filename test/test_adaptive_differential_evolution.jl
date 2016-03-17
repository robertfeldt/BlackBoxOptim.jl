NumTestRepetitions = 100

facts("Adaptive differential evolution optimizer") do

ss = symmetric_search_space(1, (0.0, 10.0))
fake_problem = FunctionBasedProblem(x -> 0.0, "test_problem", MinimizingFitnessScheme, ss)

ade = adaptive_de_rand_1_bin(fake_problem, ParamsDict(
  :Population => rand(1, 100)))

context("parameters adjust!()") do
  for(i in 1:NumTestRepetitions)
    cur_cr, cur_f = BlackBoxOptim.crossover_parameters(ade.modify, ade.population, rand(1:popsize(ade.population)))
    @fact (0.0 <= cur_cr <= 1.0) --> true
    @fact (0.0 <= cur_f <= 1.0) --> true
    BlackBoxOptim.adjust!(ade.modify, 0, 1, 0.0, 0.0, false)
    # FIXME this fails too often; needs adjusting the distribution?
    #new_cr, new_f = BlackBoxOptim.crossover_parameters(ade.params, 1)
    #@fact new_cr != cur_cr --> true
    #@fact new_f != cur_f --> true
  end
end

# FIXME actually this test is not required as the standard DE already tests for that
context("ask()") do
  for(i in 1:NumTestRepetitions)
    res = BlackBoxOptim.ask(ade)
    @fact length(res) --> 2

    trial, target = res

    @fact ndims(trial.params) --> 1
    @fact (1 <= trial.index <= popsize(ade)) --> true
    @fact in(trial.params, ade.embed.searchSpace) --> true

    @fact ndims(target.params) --> 1
    @fact (1 <= target.index <= popsize(ade)) --> true
    @fact in(target.params, ade.embed.searchSpace) --> true

    @fact trial.index --> target.index
  end
end

end
