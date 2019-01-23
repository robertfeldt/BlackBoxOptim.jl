using BlackBoxOptim

function fitness_for_opt(problem, numDimensions, populationSize, numSteps, method)
    problem = instantiate(problem, numDimensions)

    println("\n$(problem.name), n = $(numdims(problem)), optimizer = $(string(method))")

    res = bboptimize(problem; Method = method,
        NumDimensions = numDimensions,
        PopulationSize = populationSize,
        MaxSteps = numSteps
    )

    return best_fitness(res)
end
