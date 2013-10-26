using GlobalOptim
using GlobalOptim.Problems

function de_fitness_for(problem, numDimensions, populationSize, numSteps)
  problem = GlobalOptim.Problems.set_numdims!(numDimensions, problem)

  ss = search_space(problem)
  pop = rand_population(populationSize, ss)
  de = DEOpt(
    pop, ss,
    {"f" => 0.4, "cr" => 0.5, "NumParents" => 3}
  )

  best, fitness = optimize(problem, de, numSteps)
  fitness
end

facts("Optimize single objective problems in 5, 10, and 30 dimensions with DE") do
  simple_problems = ["Sphere", "Schwefel2.22", "Schwefel2.22"]
  for(problem in simple_problems)
    p = GlobalOptim.Problems.examples[problem]
    @fact de_fitness_for(p, 5, 20,  2e3) < 0.01 => true
    @fact de_fitness_for(p, 10, 20, 5e3) < 0.01 => true
    @fact de_fitness_for(p, 30, 25, 2e4) < 0.01 => true
  end

  harder_problems = ["Schwefel1.2"]
  for(problem in harder_problems)
    p = GlobalOptim.Problems.examples[problem]
    @fact de_fitness_for(p, 5, 20,  5e3) < 0.01 => true
    @fact de_fitness_for(p, 10, 20, 2e4) < 0.01 => true
    # Have to investigate why DE so bad for this one:
    @fact de_fitness_for(p, 30, 20, 1e5) < 100.0 => true
  end
end
