using BlackBoxOptim

function fitness_for_opt(problem, numDimensions, populationSize, numSteps, method)

  problem = fixed_dim_problem(problem, numDimensions)

  println("\n$(problem.name), n = $(numdims(problem)), optimizer = $(string(method))")

  best, fitness = bboptimize(problem; method = method, parameters = {
    :NumDimensions => numDimensions,
    :PopulationSize => populationSize,
    :MaxSteps => numSteps
    })

  fitness
end
