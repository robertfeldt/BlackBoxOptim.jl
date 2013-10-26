abstract Population
abstract PopulationWithFitness <: Population

type FloatVectorPopulation <: PopulationWithFitness
  # The population is a matrix of floats, each row being an individual.
  individuals::Array{Float64, 2}

  # The fitnesses is a matrix of floats, each row being the fitness for one
  # individual, and each column corresponding to one objective/goal.
  fitness::Array{Float64, 2}
  fitness_scheme::FitnessScheme

  # A population always saves the top/best individuals seen during a search,
  # together with their fitness values.
  top::Array{Float64, 2}
  top_fitness::Array{Float64, 2}
  top_size::Int

  function FloatVectorPopulation(size = 100, dimensions = 1, 
    scheme = float_vector_scheme_min(), numTopIndividuals = 10)
    inds = rand(size, dimensions)
    fs = ones(size, dimensions) * scheme.worst_fitness # Bug! Need not be the same num of objectives as there are dimensions...
    new(inds, fs, scheme, inds[1:numTopIndividuals,:], fs[1:numTopIndividuals,:], numTopIndividuals)
  end
end

best(pop::FloatVectorPopulation) = pop.top[1]
bestfitness(pop::FloatVectorPopulation) = pop.top_fitness[1]

# Update the toplist if this individual is good enough.
function update_toplist!(candidate, candidateFitness, population::PopulationWithFitness)
  fscheme = population.fitness_scheme
  if isbetter(candidateFitness, top_fitness[end], scheme)
    i = length(population.top)
    fs = population.top_fitness
    while(i >= 1 & isbetter(candidateFitness, fs[i], fscheme))
      i -= 1
    end
    population.top = vcat(population.top[1:i,:], candidate, population.top[i:end-1,:])
    population.top_fitness = vcat(population.top_fitness[1:i,:], candidateFitness, population.top_fitness[i:end-1,:])
  end
end

# Latin hypercube sampling of numIndividuals in a given search space.
function latin_hypercube_sampling(numSamples, searchSpace::Vector{(Float64, Float64)})
  numdims = length(searchSpace)
  result = zeros(numSamples, numdims)
  for(i in 1:numdims)
    min, max = searchSpace[i]
    samples_for_dim = linspace(min, max, numSamples)
    result[:,i] = samples_for_dim[shuffle([1:numSamples])]
  end
  result
end