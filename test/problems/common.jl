using BlackBoxOptim
using BlackBoxOptim.Problems

function fitness_for_opt(problem, numDimensions, populationSize, numSteps, 
  optFunc = de_rand_1_bin_radiuslimited)
  problem = BlackBoxOptim.Problems.set_numdims!(numDimensions, problem)

  ss = search_space(problem)

  pop = BlackBoxOptim.rand_individuals_lhs(ss, populationSize)

  opt = optFunc(ss; population = pop)

  println("\n$(problem.name), n = $(numdims(problem)), optimizer = $(name(opt))")

  best, fitness = BlackBoxOptim.run_optimizer_on_problem(opt, problem, numSteps)
  fitness
end
