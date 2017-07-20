# The paper:
#
#  Yi Sun, F. Gomez, T. Schaul, J. Schmidhuber. A Linear Time Natural Evolution
#  Strategy for Non-Separable Functions. In Proceedings of the Genetic and
#  Evolutionary Computation Conference (GECCO), Amsterdam, 2013.
#
# proposes a variant of a Natural Evolution Strategy which can successfully
# solve very high-dimensional problems. In particular they show state-of-the-art
# performance on the 512-dimensional Rosenbrock function.
#
# Here we try our DE implementations on this problem to see what kind of
# performance we get.

using BlackBoxOptim

function fitness_for_opt(problemname, numDimensions, populationSize, numSteps, method)

  p = BlackBoxOptim.example_problems[problemname]
  problem = instantiate(p, numDimensions)

  println("\n$(problem.name), n = $(numdims(problem)), optimizer = $(string(method))")

  result = bboptimize(problem; Method = method, 
    NumDimensions = numDimensions, PopulationSize = populationSize,
    MaxSteps = numSteps
  )

  best_fitness(result)
end

problemname = "Rosenbrock"

# Call once to compile everything
fitness_for_opt(problemname, 2, 10, 100, :adaptive_de_rand_1_bin)

for (dims, steps) in [(128, 1e6), (256, 2e6), (512, 4e6)]
  println("Num dims = ", dims, ", Num steps = ", steps)
  fitness = fitness_for_opt(problemname, dims, 50, steps, :adaptive_de_rand_1_bin)
end
