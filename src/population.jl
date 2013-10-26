abstract Population
abstract PopulationWithFitness <: Population

type FloatVectorPopulation <: PopulationWithFitness
  # The population is a matrix of floats, each row being an individual.
  inidividuals::Array{Float64, 2}

  # The fitnesses is a matrix of floats, each row being the fitness for one
  # individual, and each column corresponding to one objective/goal.
  fitness::Array{Float64, 2}

  # The currently best individual, and its fitness, is also saved here.
  best::Array{Float64, 2}
  best_fitness::Array{Float64, 2}
  best_index::Int
end

# Compare fitness by the min of their sum, lower is better.
function minimize_fitness_sum_comparator(f1::Array{Float64, 2}, f2::Array{Float64, 2})
  sum(f1) <= sum(f2)
end

# Update the best individual if this one is better.
function update_best(candidate, candidateFitness, population::PopulationWithFitness, 
  better_than_comparator = minimize_fitness_sum_comparator)
  if(better_than_comparator(candidateFitness, population.best_fitness))
  end
end