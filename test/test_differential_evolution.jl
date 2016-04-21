facts("Differential evolution optimizer") do

ss = symmetric_search_space(1, (0.0, 10.0))
fake_problem = FunctionBasedProblem(x -> 0.0, "test_problem", MinimizingFitnessScheme, ss)
DE = de_rand_1_bin(fake_problem, ParamsDict(
  :Population => collect(1.0:10.0)',
  :f => 0.4, :cr => 0.5, :NumParents => 3))

context("ask()/tell!()") do
  for i in 1:NumTestRepetitions
    res = BlackBoxOptim.ask(DE)
    @fact BlackBoxOptim.candi_pool_size(BlackBoxOptim.population(DE)) --> 0 # no candidates in the pool as we just exhausted it
    @fact length(res) --> 2

    trial, target = res

    @fact ndims(trial.params) --> 1
    @fact (1 <= trial.index <= popsize(DE.population)) --> true
    @fact in(trial.params, DE.embed.searchSpace) --> true

    @fact ndims(target.params) --> 1
    @fact (1 <= target.index <= popsize(DE.population)) --> true
    @fact in(target.params, DE.embed.searchSpace) --> true

    @fact trial.index --> target.index

    BlackBoxOptim.tell!(DE, res)
    @fact BlackBoxOptim.candi_pool_size(BlackBoxOptim.population(DE)) --> 2 # test that all candidates returned to the pool
  end
end

end
