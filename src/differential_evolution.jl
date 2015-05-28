using Distributions

abstract DifferentialEvolutionOpt <: PopulationOptimizer

name(de::DifferentialEvolutionOpt) = de.name

DE_DefaultOptions = {
  "f" => 0.6,
  "cr" => 0.7,
  "NumParents" => 3,
  "SamplerRadius" => 8,
}

type DiffEvoOpt <: DifferentialEvolutionOpt
  name::ASCIIString

  # A population is a matrix of floats.
  population::Array{Float64, 2}

  # A search space is defined by the min and max values (in tuples) for each
  # of its dimenions. The dimension is the length of an individual, i.e. the
  # number of Float64 values in it.
  search_space::SearchSpace

  # Options
  options

  # Set of functions that together define a specific DE strategy.
  sample::Function
  mutate::Function
  crossover::Function
  bound::Function

  function DiffEvoOpt(name, pop, ss, options, sample, mutate, crossover, bound)
    new(name, pop, ss, mergeparam(DE_DefaultOptions, options), sample, mutate, crossover, bound)
  end
end

popsize(opt::DifferentialEvolutionOpt) = Base.size(population(opt),1)
fconst(de::DifferentialEvolutionOpt, i) = de.options["f"]
crconst(de::DifferentialEvolutionOpt, i) = de.options["cr"]

# Ask for a new candidate object to be evaluated, and a list of individuals
# it should be ranked with. The individuals are supplied as an array of tuples
# with the individual and its index.
function ask(de::DifferentialEvolutionOpt)
  # Sample parents and target
  numparents = de.options["NumParents"]
  indices = de.sample(de, 1 + numparents)
  #print("indices = "); show(indices); println("")
  parent_indices = indices[1:numparents]
  #print("parent_indices = "); show(parent_indices); println("")
  target_index = indices[end]
  #print("target_index = "); show(target_index); println("")
  target = de.population[target_index,:]
  #print("target = "); show(target); println("")

  donor = de.mutate(de, target_index, parent_indices)
  #print("donor = "); show(donor); println("")

  # Crossover donor and target
  trial = de.crossover(de, target, target_index, donor)

  # Bound the trial vector according to search space bounds
  trial = de.bound(trial, target, de.search_space)

  # Return the candidates that should be ranked as tuples including their
  # population indices.
  return [(trial, target_index), (target, target_index)]
end

function random_sampler(de::DifferentialEvolutionOpt, numSamples)
  sample(1:popsize(de), numSamples; replace = false)
end

# This implements a "trivial geography" similar to Spector and Kline (2006)
# by first sampling an individual randomly and then selecting additional
# individuals for the same tournament within a certain deme of limited size
# for the sub-sequent individuals in the population. The version we implement
# here is from:
#  I. Harvey, "The Microbial Genetic Algorithm", in Advances in Artificial Life
#  Darwin Meets von Neumann, Springer, 2011.
# The original paper is:
#  Spector, L., and J. Klein. 2005. Trivial Geography in Genetic Programming.
#  In Genetic Programming Theory and Practice III, edited by T. Yu, R.L. Riolo,
#  and B. Worzel, pp. 109-124. Boston, MA: Kluwer Academic Publishers.
#  http://faculty.hampshire.edu/lspector/pubs/trivial-geography-toappear.pdf
#
function radius_limited_sampler(de::DifferentialEvolutionOpt, numSamples)
  # The radius must be at least as big as the number of samples + 2 so that
  # there is something to sample from.
  radius = max(de.options["SamplerRadius"], numSamples+2)
  psize = popsize(de)
  deme_start = rand(1:psize)
  indices = sample(deme_start:(deme_start+radius-1), numSamples; replace = false)
  # Ensure they are not out of bounds by wrapping over at the end.
  map(indices) do index
    if index > psize
      mod(index, psize) + 1 # We have to increase by 1 since Julia arrays start indices at 1
    else
      index
    end
  end
end

# DE/rand/1 mutation strategy
function de_mutation_rand_1(de::DifferentialEvolutionOpt, targetIndex, parentIndices)
  f = fconst(de, targetIndex)
  p = de.population[parentIndices,:]
  return p[3,:] + (f * (p[1,:] - p[2,:]))
end

# DE/rand/2 mutation strategy
function de_mutation_rand_2(de::DifferentialEvolutionOpt, targetIndex, parentIndices)
  f = fconst(de, targetIndex)
  p = de.population[parentIndices,:]
  return p[3,:] + (f * (p[1,:] - p[2,:])) + (f * (p[4,:] - p[5,:]))
end

# Binomial crossover for DE, i.e. DE/*/*/bin.
function de_crossover_binomial(de::DifferentialEvolutionOpt, target, targetIndex, donor)
  trial = copy(target)

  # Always ensure at least one value from donor is copied to trial vector
  jrand = rand(1:length(trial))
  trial[jrand] = donor[jrand]

  # Now crossover randomly for the rest of the indices
  switch = rand(length(trial)) .<= crconst(de, targetIndex)
  #print("switch = "); show(switch); println("")
  #print("trial = "); show(trial); println("")
  #print("donor = "); show(donor); println("")
  trial[:,switch] = donor[:,switch]
  #print("trial = "); show(trial); println("")

  return trial
end

# If we come out-of-bounds we randomly sample between the target value
# and the bound.
function rand_bound_from_target!(individual, target, searchSpace::SearchSpace)
  ssmins, ssmaxs = mins(searchSpace), maxs(searchSpace)
  for i in 1:length(individual)
    min, max = ssmins[i], ssmaxs[i]

    if individual[i] < min
      individual[i] = min + rand() * (target[i] - min)
    elseif individual[i] > ssmaxs[i]
      individual[i] = target[i] + rand() * (max - target[i])
    end
  end
  individual
end

function rand_bound_from_target!(individual, target, searchSpace::Array{(Float64,Float64),1})
  rand_bound_from_target!(individual, target, RangePerDimSearchSpace(searchSpace))
end

# Tell the optimizer about the ranking of candidates. Returns the number of
# better candidates that were inserted into the population.
function tell!(de::DiffEvoOpt,
  # archive::Archive, # Skip for now
  rankedCandidates)
  num_candidates = length(rankedCandidates)
  num_better = 0
  for i in 1:div(num_candidates, 2)
    candidate, index = rankedCandidates[i]
    if candidate != de.population[index, :]
      num_better += 1
      #print("candidate = "); show(candidate); println("")
      #print("index = "); show(index); println("")
      #print("target = "); show(de.population[index,:]); println("")
      old = de.population[index,:]
      de.population[index,:] = candidate
    end
  end
  num_better
end

# Now we can create specific DE optimizers that are commonly used in the
# literature.

# The most used DE/rand/1/bin.
function de_rand_1_bin(parameters = Dict(); sampler = random_sampler, name = "DE/rand/1/bin")
  params = Parameters(parameters, DE_DefaultOptions)
  # Ensure NumParents is 3 since de_mutation_rand_1 requires it.
  params["NumParents"] = 3
  DiffEvoOpt(name, params[:Population], params[:SearchSpace], params,
    sampler,
    de_mutation_rand_1,
    de_crossover_binomial,
    rand_bound_from_target!)
end

function de_rand_2_bin(parameters = Dict(); sampler = random_sampler, name = "DE/rand/2/bin")
  params = Parameters(parameters, DE_DefaultOptions)
  # Ensure NumParents is 5 since de_mutation_rand_2 requires it.
  params["NumParents"] = 5
  DiffEvoOpt(name, params[:Population], params[:SearchSpace], params,
    sampler,
    de_mutation_rand_2,
    de_crossover_binomial,
    rand_bound_from_target!)
end

# The most used DE/rand/1/bin with "local geography" via radius limited sampling.
function de_rand_1_bin_radiuslimited(parameters = Dict())
  de_rand_1_bin(parameters; sampler = radius_limited_sampler, name = "DE/rand/1/bin/radiuslimited")
end

function de_rand_2_bin_radiuslimited(parameters = Dict())
  de_rand_2_bin(parameters; sampler = radius_limited_sampler, name = "DE/rand/2/bin/radiuslimited")
end
