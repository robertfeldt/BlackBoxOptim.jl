DE = de_rand_1_bin(@compat Dict{Symbol,Any}(
  :SearchSpace => symmetric_search_space(1, (0.0, 10.0)),
  :Population => reshape(collect(1.0:10.0), 10, 1),
  :f => 0.4,
  :cr => 0.5,
  :NumParents => 3))

facts("Differential evolution optimizer") do

context("random_sampler") do
  for(i in 1:NumTestRepetitions)
    numSamples = rand(1:8)
    sampled = BlackBoxOptim.random_sampler(DE, numSamples)

    @fact length(sampled) => numSamples

    popindices = 1:popsize(DE)

    # All sampled indices are indices into the population
    @fact all([in(index, popindices) for index in sampled]) => true
  end
end

context("radius_limited_sampler") do
  DE = de_rand_1_bin(@compat Dict{Symbol,Any}(
    :SearchSpace => symmetric_search_space(1, (0.0, 10.0)),
    :Population => rand(100,1),
    :f => 0.4, :cr => 0.5, :NumParents => 3))

  psize = popsize(DE)
  popindices = 1:popsize(DE)
  sampler_radius = DE.options[:SamplerRadius]

  for(i in 1:NumTestRepetitions)
    numSamples = rand(1:sampler_radius)
    sampled = BlackBoxOptim.radius_limited_sampler(DE, numSamples)

    @fact length(sampled) => numSamples

    # All sampled indices are indices into the population
    @fact all([in(index, popindices) for index in sampled]) => true

    mini, maxi = minimum(sampled), maximum(sampled)
    if (maxi - mini) > max(numSamples+2, sampler_radius)
      @fact mini + psize <= maxi + sampler_radius => true
    end
  end
end

context("rand_bound_from_target!") do
  context("does nothing if within bounds") do
    @fact BlackBoxOptim.rand_bound_from_target!([0.0], [0.0], [(0.0, 1.0)]) => [0.0]

    @fact BlackBoxOptim.rand_bound_from_target!([0.0, 11.4], [0.1, 12.3],
      RangePerDimSearchSpace([(0.0, 1.0), (10.0, 15.0)])) => [0.0, 11.4]
  end

  context("bounds if lower than min bound") do
    @fact BlackBoxOptim.rand_bound_from_target!([-0.1], [0.0], [(0.0, 1.0)]) => [0.0]

    res = BlackBoxOptim.rand_bound_from_target!([-0.1], [0.5], [(0.0, 1.0)])
    @fact 0.5 >= res[1] >= 0.0 => true

    res = BlackBoxOptim.rand_bound_from_target!([-11.1, 0.5], [-10.8, 0.5], [(-11.0, 1.0), (0.0, 1.0)])
    @fact -10.8 >= res[1] >= -11.0 => true
    @fact res[2] => 0.5

    res = BlackBoxOptim.rand_bound_from_target!([50.4, -103.1], [49.6, -101.4], [(30.0, 60.0), (-102.0, -1.0)])
    @fact res[1] => 50.4
    @fact -101.4 >= res[2] >= -102.0 => true
  end

  context("bounds if higher than max bound") do
    @fact BlackBoxOptim.rand_bound_from_target!([1.1], [1.0], [(0.0, 1.0)]) => [1.0]

    res = BlackBoxOptim.rand_bound_from_target!([97.0], [95.0], [(-10.0, 96.0)])
    @fact 95.0 <= res[1] <= 96.0 => true
  end
end

context("de_crossover_binomial") do
  context("always copies from donor if length is 1") do
    @fact BlackBoxOptim.de_crossover_binomial(DE, [0.0], 1, [1.0]) => [1.0]

    @fact BlackBoxOptim.de_crossover_binomial(DE, [-10.0], 1, [42.42]) => [42.42]
  end

  context("always copies at least one element from donor") do
    for(i in 1:NumTestRepetitions)
      len = rand(1:100)
      target, donor = rand(1, len), rand(1, len)
      res = BlackBoxOptim.de_crossover_binomial(DE, target, 1, donor)
      @fact any([ in(x, donor) for x = res ]) => true
    end
  end

  context("unlikely to copy everything if vectors are large") do
    for(i in 1:NumTestRepetitions)
      len = 50
      target, donor = rand(1, len), rand(1, len)
      res = BlackBoxOptim.de_crossover_binomial(DE, target, 1, donor)
      @fact any([ in(x, target) for x = res ]) => true
    end
  end
end

context("de_mutation_rand_1") do
  @fact ndims(BlackBoxOptim.de_mutation_rand_1(DE, 1, [4,9,8])) => 2

  @fact BlackBoxOptim.de_mutation_rand_1(DE, 1, [1, 2, 3])[1] => (3.0 + (0.4 * (1.0 - 2.0)))
  @fact BlackBoxOptim.de_mutation_rand_1(DE, 2, [2, 3, 1])[1] => (1.0 + (0.4 * (2.0 - 3.0)))
  @fact BlackBoxOptim.de_mutation_rand_1(DE, 3, [3, 2, 1])[1] => (1.0 + (0.4 * (3.0 - 2.0)))
  @fact BlackBoxOptim.de_mutation_rand_1(DE, 4, [1, 3, 5])[1] => (5.0 + (0.4 * (1.0 - 3.0)))
  @fact BlackBoxOptim.de_mutation_rand_1(DE, 5, [4, 9, 8])[1] => (8.0 + (0.4 * (4.0 - 9.0)))

  de2 = de_rand_1_bin(@compat Dict{Symbol,Any}(:SearchSpace => symmetric_search_space(2, (0.0, 10.0)),
    :Population => reshape(collect(1.0:8.0), 4, 2),
    :f => 0.6, :cr => 0.5, :NumParents => 3)
  )

  res = BlackBoxOptim.de_mutation_rand_1(de2, 10, [1,2,3])
  @fact ndims(res) => 2
  @fact res[1] => (3.0 + (0.6 * (1.0 - 2.0)))
  @fact res[2] => (7.0 + (0.6 * (5.0 - 6.0)))

  res2 = BlackBoxOptim.de_mutation_rand_1(de2, 2, [1,2,4])
  @fact res2[1] => (4.0 + (0.6 * (1.0 - 2.0)))
  @fact res2[2] => (8.0 + (0.6 * (5.0 - 6.0)))
end

context("ask") do
  for(i in 1:NumTestRepetitions)
    res = BlackBoxOptim.ask(DE)

    @fact length(res) => 2
    trial, trialIndex = res[1]
    target, targetIndex = res[2]

    @fact ndims(trial) => 2
    @fact 1 <= trialIndex <= length(DE.population) => true
    @fact isinspace(trial, DE.search_space) => true

    @fact ndims(target) => 2
    @fact 1 <= targetIndex <= length(DE.population) => true
    @fact isinspace(target, DE.search_space) => true

    @fact trialIndex == targetIndex => true
  end
end

end
