NumTestRepetitions = 100

facts("rand_bound_from_target!") do
  context("does nothing if within bounds") do
    @fact GlobalOptim.rand_bound_from_target!([0.0], [0.0], [(0.0, 1.0)]) => [0.0]

    @fact GlobalOptim.rand_bound_from_target!([0.0, 11.4], [0.1, 12.3], [(0.0, 1.0), (10.0, 15.0)]) => [0.0, 11.4]
  end

  context("bounds if lower than min bound") do
    @fact GlobalOptim.rand_bound_from_target!([-0.1], [0.0], [(0.0, 1.0)]) => [0.0]

    res = GlobalOptim.rand_bound_from_target!([-0.1], [0.5], [(0.0, 1.0)])
    @fact 0.5 >= res[1] >= 0.0 => true

    res = GlobalOptim.rand_bound_from_target!([-11.1, 0.5], [-10.8, 0.5], [(-11.0, 1.0), (0.0, 1.0)])
    @fact -10.8 >= res[1] >= -11.0 => true
    @fact res[2] => 0.5

    res = GlobalOptim.rand_bound_from_target!([50.4, -103.1], [49.6, -101.4], [(30.0, 60.0), (-102.0, -1.0)])
    @fact res[1] => 50.4
    @fact -101.4 >= res[2] >= -102.0 => true
  end

  context("bounds if higher than max bound") do
    @fact GlobalOptim.rand_bound_from_target!([1.1], [1.0], [(0.0, 1.0)]) => [1.0]

    res = GlobalOptim.rand_bound_from_target!([97.0], [95.0], [(-10.0, 96.0)])
    @fact 95.0 <= res[1] <= 96.0 => true
  end
end

DE = DEOpt(
  reshape([1.0:10.0], 10, 1),
  [(0.0, 10.0)],
  {"f" => 0.4, "cr" => 0.5, "NumParents" => 3}
)

facts("de_crossover_binomial") do
#  context("always copies from donor if length is 1") do
#    @fact GlobalOptim.de_crossover_binomial(DE, [0.0], [1.0]) => [1.0]
#
#    @fact GlobalOptim.de_crossover_binomial(DE, [-10.0], [42.42]) => [42.42]
#  end

  context("always copies at least one element from donor") do
    for(i in 1:NumTestRepetitions)
      len = rand(1:100)
      target, donor = rand(len, 1), rand(len, 1)
      res = GlobalOptim.de_crossover_binomial(DE, target, donor)
      @fact any([ in(x, donor) for x = res ]) => true
    end
  end

  context("unlikely to copy everything if vectors are large") do
    for(i in 1:NumTestRepetitions)
      len = 50
      target, donor = rand(len, 1), rand(len, 1)
      res = GlobalOptim.de_crossover_binomial(DE, target, donor)
      @fact any([ in(x, target) for x = res ]) => true
    end
  end
end

facts("de_mutation_rand_1") do
  @fact ndims(GlobalOptim.de_mutation_rand_1(DE, [4,9,8])) => 2

  @fact GlobalOptim.de_mutation_rand_1(DE, [1, 2, 3])[1] => (3.0 + (0.4 * (1.0 - 2.0)))
  @fact GlobalOptim.de_mutation_rand_1(DE, [2, 3, 1])[1] => (1.0 + (0.4 * (2.0 - 3.0)))
  @fact GlobalOptim.de_mutation_rand_1(DE, [3, 2, 1])[1] => (1.0 + (0.4 * (3.0 - 2.0)))
  @fact GlobalOptim.de_mutation_rand_1(DE, [1, 3, 5])[1] => (5.0 + (0.4 * (1.0 - 3.0)))
  @fact GlobalOptim.de_mutation_rand_1(DE, [4, 9, 8])[1] => (8.0 + (0.4 * (4.0 - 9.0)))

  DE2 = DEOpt(
    reshape([1.0:8.0], 4, 2),
    [(0.0, 10.0), (0.0, 10.0)],
    {"f" => 0.4, "cr" => 0.5, "NumParents" => 3}
  )
  res = GlobalOptim.de_mutation_rand_1(DE2, [1,2,3])
  @fact ndims(res) => 2
  @fact res[1] => (3.0 + (0.4 * (1.0 - 2.0)))
  @fact res[2] => (7.0 + (0.4 * (5.0 - 6.0)))

  res2 = GlobalOptim.de_mutation_rand_1(DE2, [1,2,4])
  @fact res2[1] => (4.0 + (0.4 * (1.0 - 2.0)))
  @fact res2[2] => (8.0 + (0.4 * (5.0 - 6.0)))
end

function isWithinBounds(ind, de)
  ss = de.search_space
  all([(ss[i][1] <= ind[i] <= ss[i][2]) for i=1:length(ind)])
end

facts("ask") do
  for(i in 1:NumTestRepetitions)
    res = GlobalOptim.ask(DE)

    @fact length(res) => 2
    trial, trialIndex = res[1]
    target, targetIndex = res[2]

    @fact ndims(trial) => 2
    @fact 1 <= trialIndex <= length(DE.population) => true
    @fact isWithinBounds(trial, DE) => true

    @fact ndims(target) => 2
    @fact 1 <= targetIndex <= length(DE.population) => true
    @fact isWithinBounds(target, DE) => true

    @fact trialIndex == targetIndex => true
  end
end