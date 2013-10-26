using GlobalOptim
using GlobalOptim.Problems

function de_fitness_for(problem, numDimensions, populationSize, numSteps)
  problem = GlobalOptim.Problems.set_numdims!(numDimensions, problem)
  println("\n$(problem.name), n = $(numdims(problem))")

  ss = search_space(problem)
  pop = rand_population(populationSize, ss)
  de = de_rand_1_bin(pop, ss)

  best, fitness = optimize(problem, de, numSteps)
  fitness
end

facts("Optimize single objective problems in 5, 10, and 30 dimensions with DE") do
  simple_problems = ["Sphere", "Schwefel2.22", "Schwefel2.22"]
  for(problem in simple_problems)
    context(problem) do
      p = GlobalOptim.Problems.examples[problem]
      @fact de_fitness_for(p, 5, 20,  10e3) < 0.01 => true
      @fact de_fitness_for(p, 10, 20, 5e4) < 0.01 => true
      @fact de_fitness_for(p, 30, 25, 1e5) < 0.01 => true
    end
  end

  context("Schwefel1.2") do
    problem = "Schwefel1.2"
    p = GlobalOptim.Problems.examples[problem]
    @fact de_fitness_for(p, 5, 20,  5e3) < 0.01 => true
    @fact de_fitness_for(p, 10, 50, 1e5) < 10.0 => true
    # Have to investigate why DE so bad for this one:
    #@fact de_fitness_for(p, 30, 100, 3e5) < 100.0 => true
  end

  context("Rosenbrock") do
    problem = "Rosenbrock"
    p = GlobalOptim.Problems.examples[problem]
    @fact de_fitness_for(p, 5, 20,   1e4) < 100.0 => true
    @fact de_fitness_for(p, 10, 20,  5e4) < 100.0 => true
    @fact de_fitness_for(p, 30, 40, 1e5) < 100.0 => true
  end
end