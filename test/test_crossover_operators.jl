facts("Crossover operators") do

ss = symmetric_search_space(1, (0.0, 10.0))
fake_problem = FunctionBasedProblem(x -> 0.0, "test_problem", MinimizingFitnessScheme, ss)
DE = de_rand_1_bin(fake_problem, ParamsDict(
  :Population => collect(1.0:10.0)',
  :f => 0.4, :cr => 0.5, :NumParents => 3))

context("DiffEvoRandBin1") do
  context("always copies from donor if length is 1") do
    @fact BlackBoxOptim.apply!(BlackBoxOptim.DiffEvoRandBin1(0.0, 0.0),
                               [0.0], 4, DE.population, [1,2,3]) --> [3.0]
  end

  context("always copies at least one element from donor") do
    for i in 1:NumTestRepetitions
      len = rand(1:100)
      pop = rand(len,4)
      target = pop[:,1]
      saved_target = copy(target)
      res = BlackBoxOptim.apply!(BlackBoxOptim.DiffEvoRandBin1(0.0, 0.0),
                                target, 1, pop, [2,3,4])
      @fact (res === target) --> true
      @fact ndims( res ) --> 1
      @fact length( res ) --> len
      @fact sum(target .!= saved_target) --> 1
    end
  end

  context("unlikely to copy everything if vectors are large") do
    for i in 1:NumTestRepetitions
      len = 50 + rand(1:50)
      pop = rand(len,4)
      target = pop[:,1]
      saved_target = copy(target)
      res = BlackBoxOptim.apply!(BlackBoxOptim.DiffEvoRandBin1(0.1, 0.5),
                                 target, 1, pop, [2,3,4])
      @fact any(target .!= saved_target) --> true
      @fact any(target .== saved_target) --> true
    end
  end

  context("correctly modifies the parameters vector") do
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                [0.0], 4, DE.population, [1,2,3]) --> [3.0 + (0.4 * (1.0 - 2.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                [0.0], 4, DE.population, [2,3,1]) --> [1.0 + (0.4 * (2.0 - 3.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                [0.0], 4, DE.population, [3,2,1]) --> [1.0 + (0.4 * (3.0 - 2.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                [0.0], 4, DE.population, [1,3,5]) --> [5.0 + (0.4 * (1.0 - 3.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.4),
                                [0.0], 3, DE.population, [4,9,8]) --> [8.0 + (0.4 * (4.0 - 9.0))]

    pop2 = reshape(collect(1.0:8.0), 4, 2)'
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.6),
                                [0.0,0.0], 4, pop2, [1,2,3]) --> [3.0 + (0.6 * (1.0 - 2.0)),
                                                                 7.0 + (0.6 * (5.0 - 6.0))]
    @fact BlackBoxOptim.apply!( BlackBoxOptim.DiffEvoRandBin1(1.0, 0.6),
                                [0.0,0.0], 4, pop2, [1,2,4]) --> [4.0 + (0.6 * (1.0 - 2.0)),
                                                                 8.0 + (0.6 * (5.0 - 6.0))]
  end
end

end
