include("bimodal_cauchy_distribution.jl")

ADE_DefaultOptions = merge(DE_DefaultOptions, {
  # Distributions we will use to generate new F and CR values.
  "fdistr" => bimodal_cauchy(0.65, 0.1, 1.0, 0.1),
  "crdistr" => bimodal_cauchy(0.1, 0.1, 0.95, 0.1),
})

# An Adaptive DE typically change parameters of the search dynamically. This is 
# typically done in the tell! function when we know if the trial vector
# was better than the target vector.
type AdaptConstantsDiffEvoOpt <: DifferentialEvolutionOpt
  name::ASCIIString

  # A population is a matrix of floats.
  population::Array{Float64, 2}

  search_space::SearchSpace

  # Options
  options::Dict{Any,Any}

  # Set of functions that together define a specific DE strategy.
  sample::Function
  mutate::Function
  crossover::Function
  bound::Function

  # Specific data and functions for adaptation
  fs::Vector{Float64}   # One f value per individual in population
  crs::Vector{Float64}  # One cr value per individual in population

  function AdaptConstantsDiffEvoOpt(name, pop, ss, options, sample, mutate, crossover, bound)
    popsize = size(pop, 1)
    fs = [sample_bimodal_cauchy(options["fdistr"]; truncateBelow0 = false) for i in 1:popsize]
    crs = [sample_bimodal_cauchy(options["crdistr"]) for i in 1:popsize]
    new(name, pop, ss, merge(DE_DefaultOptions, options), 
      sample, mutate, crossover, bound, fs, crs)
  end
end

# To get the constants for an adaptive DE we access the vectors of constants.
fconst(ade::AdaptConstantsDiffEvoOpt, i) = ade.fs[i]
crconst(ade::AdaptConstantsDiffEvoOpt, i) = ade.crs[i]

# To sample we use the distribution given as options
sample_f(ade::AdaptConstantsDiffEvoOpt) = sample_bimodal_cauchy(ade.options["fdistr"]; truncateBelow0 = false)
sample_cr(ade::AdaptConstantsDiffEvoOpt) = sample_bimodal_cauchy(ade.options["crdistr"]; truncateBelow0 = false)

# Tell the optimizer about the ranking of candidates. Returns the number of
# better candidates that were inserted into the population.
function tell!(de::AdaptConstantsDiffEvoOpt, 
  # archive::Archive, # Skip for now
  rankedCandidates)
  num_candidates = length(rankedCandidates)
  num_better = 0
  for i in 1:div(num_candidates, 2)
    candidate, index = rankedCandidates[i]
    if candidate != de.population[index, :]
      num_better += 1
      old = de.population[index,:]
      de.population[index,:] = candidate
      # Since the trial vector was better we keep the f and cr values for this target.
    else
      # The trial vector for this target was not better so we change the f and cr constants.
      de.fs[index] = sample_f(de)
      de.crs[index] = sample_cr(de)
    end
  end
  num_better
end

function adaptive_de_rand_1_bin(population = rand(100,1), 
  searchSpace = RangePerDimSearchSpace([(0.0, 1.0)]), options = ADE_DefaultOptions)
  AdaptConstantsDiffEvoOpt("AdaptiveDE/rand/1/bin", population, searchSpace, options, 
    random_sampler, 
    de_mutation_rand_1, 
    de_crossover_binomial, 
    rand_bound_from_target!)
end

function adaptive_de_rand_1_bin_radiuslimited(population = rand(100,1), 
  searchSpace = RangePerDimSearchSpace([(0.0, 1.0)]), options = ADE_DefaultOptions)
  AdaptConstantsDiffEvoOpt("AdaptiveDE/rand/1/bin/radiuslimited", population, searchSpace, options, 
    radius_limited_sampler, 
    de_mutation_rand_1, 
    de_crossover_binomial, 
    rand_bound_from_target!)
end
