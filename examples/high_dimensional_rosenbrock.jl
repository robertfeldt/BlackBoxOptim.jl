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
using BlackBoxOptim.Problems

function fitness_for_opt(problem, numDimensions, populationSize, numSteps, 
  optFunc = de_rand_1_bin_radiuslimited)
  problem = BlackBoxOptim.Problems.set_numdims!(numDimensions, problem)

  ss = search_space(problem)

  pop = BlackBoxOptim.rand_individuals_lhs(ss, populationSize)

  opt = optFunc(ss; population = pop)

  println("\n$(problem.name), n = $(numdims(problem)), optimizer = $(opt.name)")

  best, fitness = BlackBoxOptim.run_optimizer_on_problem(opt, problem, numSteps)
  fitness
end

problem = "Rosenbrock"
p = BlackBoxOptim.Problems.examples[problem]

# Call once to compile everything
fitness_for_opt(p, 2, 10, 100, adaptive_de_rand_1_bin)

for( (dims, steps) in [(128, 1e6), (256, 2e6), (512, 4e6)])
  println("Num dims = ", dims, ", Num steps = ", steps)
  fitness = fitness_for_opt(p, dims, 50, steps, adaptive_de_rand_1_bin)
end
