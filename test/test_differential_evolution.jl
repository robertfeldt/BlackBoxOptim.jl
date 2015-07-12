facts("Differential evolution optimizer") do

ss = symmetric_search_space(1, (0.0, 10.0))
fake_problem = convert(FunctionBasedProblem, x -> 0.0, "test_problem",
                       MinimizingFitnessScheme, ss) # FIXME v0.3 workaround
DE = de_rand_1_bin(fake_problem, @compat Dict{Symbol,Any}(
  :Population => collect(1.0:10.0)',
  :f => 0.4, :cr => 0.5, :NumParents => 3))

context("SimpleSelector") do
  @fact popsize(DE) => 10
  sel = BlackBoxOptim.SimpleSelector()
  for(i in 1:NumTestRepetitions)
    numSamples = rand(1:8)
    sampled = BlackBoxOptim.select(sel, DE.population, numSamples)

    @fact length(sampled) => numSamples

    # All sampled indices are indices into the population
    @fact all([in(index, 1:popsize(DE)) for index in sampled]) => true
  end
end

context("RadiusLimitedSelector") do
  local DE = de_rand_1_bin_radiuslimited(fake_problem, @compat Dict{Symbol,Any}(
    :Population => rand(1,100),
    :f => 0.4, :cr => 0.5, :NumParents => 3))

  @fact popsize(DE) => 100
  sel = DE.select
  @fact typeof(sel) => BlackBoxOptim.RadiusLimitedSelector

  for(i in 1:NumTestRepetitions)
    numSamples = rand(1:sel.radius)
    sampled = BlackBoxOptim.select(sel, DE.population, numSamples)

    @fact length(sampled) => numSamples

    # All sampled indices are indices into the population
    @fact all([in(index, 1:popsize(DE)) for index in sampled]) => true

    mini, maxi = minimum(sampled), maximum(sampled)
    if (maxi - mini) > max(numSamples+2, sel.radius)
      @fact mini + popsize(DE) <= maxi + sel.radius => true
    end
  end
end

context("RandomBound") do
  context("does nothing if within bounds") do
    @fact BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                               [0.0], [0.0]) => [0.0]

    @fact BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0), (10.0, 15.0)]),
                               [0.0, 11.4], [0.1, 12.3] ) => [0.0, 11.4]
  end

  context("bounds if lower than min bound") do
    @fact BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                               [-0.1], [0.0]) => [0.0]

    res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                               [-0.1], [0.5])
    @fact 0.5 >= res[1] >= 0.0 => true

    res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(-11.0, 1.0), (0.0, 1.0)]),
                               [-11.1, 0.5], [-10.8, 0.5])
    @fact -10.8 >= res[1] >= -11.0 => true
    @fact res[2] => 0.5

    res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(30.0, 60.0), (-102.0, -1.0)]),
                               [50.4, -103.1], [49.6, -101.4])
    @fact res[1] => 50.4
    @fact -101.4 >= res[2] >= -102.0 => true
  end

  context("bounds if higher than max bound") do
    @fact BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(0.0, 1.0)]),
                               [1.1], [1.0]) => [1.0]

    res = BlackBoxOptim.apply!(BlackBoxOptim.RandomBound([(-10.0, 96.0)]),
                               [97.0], [95.0])
    @fact 95.0 <= res[1] <= 96.0 => true
  end
end

context("DiffEvoRandBin1") do
  context("always copies from donor if length is 1") do
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 0.0, 0.0,
                                [0.0], DE.population, [1,2,3]) => [3.0]
  end

  context("always copies at least one element from donor") do
    for(i in 1:NumTestRepetitions)
      len = rand(1:100)
      pop = rand(len,4)
      target = pop[:,1]
      saved_target = copy(target)
      res = BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 0.0, 0.0,
                                  target, pop, [2,3,4])
      @fact res === target => true
      @fact ndims( res ) => 1
      @fact length( res ) => len
      @fact sum(target .!= saved_target) => 1
    end
  end

  context("unlikely to copy everything if vectors are large") do
    for(i in 1:NumTestRepetitions)
      len = 50 + rand(1:50)
      pop = rand(len,4)
      target = pop[:,1]
      saved_target = copy(target)
      res = BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 0.1, 0.5, target, pop, [2,3,4])
      @fact any(target .!= saved_target) => true
      @fact any(target .== saved_target) => true
    end
  end

  context("correctly modifies the parameters vector") do
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 1.0, 0.4,
                                [0.0], DE.population, [1,2,3]) => [3.0 + (0.4 * (1.0 - 2.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 1.0, 0.4,
                                [0.0], DE.population, [2,3,1]) => [1.0 + (0.4 * (2.0 - 3.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 1.0, 0.4,
                                [0.0], DE.population, [3,2,1]) => [1.0 + (0.4 * (3.0 - 2.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 1.0, 0.4,
                                [0.0], DE.population, [1,3,5]) => [5.0 + (0.4 * (1.0 - 3.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 1.0, 0.4,
                                [0.0], DE.population, [4,9,8]) => [8.0 + (0.4 * (4.0 - 9.0))]

    pop2 = reshape(collect(1.0:8.0), 4, 2)'
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 1.0, 0.6,
                                [0.0,0.0], pop2, [1,2,3]) => [3.0 + (0.6 * (1.0 - 2.0)),
                                                              7.0 + (0.6 * (5.0 - 6.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(), 1.0, 0.6,
                                [0.0,0.0], pop2, [1,2,4]) => [4.0 + (0.6 * (1.0 - 2.0)),
                                                              8.0 + (0.6 * (5.0 - 6.0))]
  end
end

context("ask()/tell!()") do
  for(i in 1:NumTestRepetitions)
    res = BlackBoxOptim.ask(DE)
    @fact BlackBoxOptim.candi_pool_size(BlackBoxOptim.population(DE)) => 0 # no candidates in the pool as we just exhausted it
    @fact length(res) => 2

    trial, target = res

    @fact ndims(trial.params) => 1
    @fact 1 <= trial.index <= popsize(DE.population) => true
    @fact isinspace(trial.params, DE.embed.searchSpace) => true

    @fact ndims(target.params) => 1
    @fact 1 <= target.index <= popsize(DE.population) => true
    @fact isinspace(target.params, DE.embed.searchSpace) => true

    @fact trial.index == target.index => true

    BlackBoxOptim.tell!(DE, res)
    @fact BlackBoxOptim.candi_pool_size(BlackBoxOptim.population(DE)) => 2 # test that all candidates returned to the pool
  end
end

end
