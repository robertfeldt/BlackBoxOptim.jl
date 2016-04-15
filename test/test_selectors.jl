facts("Selection operators") do

ss = symmetric_search_space(1, (0.0, 10.0))
fake_problem = FunctionBasedProblem(x -> 0.0, "test_problem", MinimizingFitnessScheme, ss)
fake_pop = collect(1.0:10.0)'

context("SimpleSelector") do
  @fact popsize(fake_pop) --> 10
  sel = BlackBoxOptim.SimpleSelector()
  for i in 1:NumTestRepetitions
    numSamples = rand(1:8)
    sampled = BlackBoxOptim.select(sel, fake_pop, numSamples)

    @fact length(sampled) --> numSamples

    # All sampled indices are indices into the population
    @fact all(index -> in(index, 1:popsize(fake_pop)), sampled) --> true
  end
end

context("RadiusLimitedSelector") do
  sel = BlackBoxOptim.RadiusLimitedSelector(20)

  for i in 1:NumTestRepetitions
    numSamples = rand(1:sel.radius)
    sampled = BlackBoxOptim.select(sel, fake_pop, numSamples)

    @fact length(sampled) --> numSamples

    # All sampled indices are indices into the population
    @fact all(index -> in(index, 1:popsize(fake_pop)), sampled) --> true

    mini, maxi = minimum(sampled), maximum(sampled)
    if (maxi - mini) > max(numSamples+2, sel.radius)
      @fact mini + popsize(fake_pop) --> less_than_or_equal(maxi + sel.radius)
    end
  end
end

end
