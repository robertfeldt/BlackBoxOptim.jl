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

facts("Optimize sphere in 5 dimensions with DE") do
  p = GlobalOptim.Problems.examples["Sphere"]
  @fact de_fitness_for(p, 5, 20, 2000) < 0.01 => true
end

facts("Optimize sphere in 10 dimensions with DE") do
  p = GlobalOptim.Problems.examples["Sphere"]
  @fact de_fitness_for(p, 10, 20, 5000) < 0.01 => true
end

facts("Optimize sphere in 30 dimensions with DE") do
  p = GlobalOptim.Problems.examples["Sphere"]
  @fact de_fitness_for(p, 30, 25, 10000) < 0.01 => true
end

facts("Optimize schwefel2.22 in 5 dimensions with DE") do
  p = GlobalOptim.Problems.examples["Schwefel2.22"]
  @fact de_fitness_for(p, 5, 20, 2000) < 0.01 => true
end

facts("Optimize schwefel2.22 in 10 dimensions with DE") do
  p = GlobalOptim.Problems.examples["Schwefel2.22"]
  @fact de_fitness_for(p, 10, 20, 5000) < 0.01 => true
end
