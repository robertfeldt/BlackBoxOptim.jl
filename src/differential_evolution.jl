using Distributions

abstract DifferentialEvolutionOpt <: PopulationOptimizer

name(de::DifferentialEvolutionOpt) = de.name

DE_DefaultOptions = @compat Dict{Symbol,Any}(
  :f => 0.6,
  :cr => 0.7,
  :SamplerRadius => 8,
)

abstract DiffEvoParameters

type FixedDiffEvoParameters <: DiffEvoParameters
  f::Float64
  cr::Float64

  FixedDiffEvoParameters(options, popsize::Int) =  new(options[:f], options[:cr])
end

crossover_parameters(params::FixedDiffEvoParameters, index) = params.cr, params.f

function adjust!( params::FixedDiffEvoParameters, index, is_improved::Bool )
    # do nothing
end

type DiffEvoOpt{P<:DiffEvoParameters,M<:MutationOperator,X<:DiffEvoCrossoverOperator,E<:EmbeddingOperator} <: DifferentialEvolutionOpt
  # TODO when sampler and bound would be parameterized, name is no longer required -- as everything is seen in the type name
  name::ASCIIString

  # A population is a matrix of floats, individuals stored in columns.
  population::Array{Float64, 2}

  # Options
  options::Parameters

  # Set of functions that together define a specific DE strategy.
  params::P        # adjust crossover parameters after fitness calculation
  sample::Function # TODO change function to SamplingOperator
  mutate::M        # mutation operator
  crossover::X     # crossover operator
  embed::E         # embedding operator
end

function DiffEvoOpt{P<:DiffEvoParameters,M<:MutationOperator,X<:DiffEvoCrossoverOperator,E<:EmbeddingOperator}(
      name::ASCIIString, pop, options::Dict{Symbol,Any}, params::P, sample,
      mutate::M = M(), crossover::X = X(), embed::E = E())
  DiffEvoOpt{P,M,X,E}(name, pop, options, params, sample, mutate, crossover, embed)
end

popsize(opt::DifferentialEvolutionOpt) = size(population(opt),2)

# Ask for a new candidate object to be evaluated, and a list of individuals
# it should be ranked with. The individuals are supplied as an array of tuples
# with the individual and its index.
function ask(de::DifferentialEvolutionOpt)
  # Sample parents and target
  indices = de.sample(de, 1 + numparents(de.crossover))
  parent_indices = indices[1:numparents(de.crossover)]
  #println("parent_indices = $(parent_indices)")
  target_index = indices[end]
  target = de.population[:,target_index] # FIXME we can cache target somewhere
  #println("target[$(target_index)] = $(target)")
  apply!( de.mutate, target )

  # Crossover parents and target
  @assert numchildren(de.crossover)==1
  trial = copy(target) # FIXME reuse some trial vector
  apply!( de.crossover,
          crossover_parameters( de.params, target_index )...,
          trial, de.population, parent_indices )
  # embed the trial parameter vector into the search space
  apply!(de.embed, trial, de.population, [target_index])
  #println("trial = $(trial)")

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
  radius = max(de.options[:SamplerRadius], numSamples+2)
  psize = popsize(de)
  deme_start = rand(1:psize)
  indices = sample(deme_start:(deme_start+radius-1), numSamples; replace = false)
  # Ensure they are not out of bounds by wrapping over at the end.
  map(indices) do index
    if index > psize
      mod(index-1, psize) + 1 # We have to increase by 1 since Julia arrays start indices at 1
    else
      index
    end
  end
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
    is_inserted = candidate != de.population[:,index]
    adjust!( de.params, index, is_inserted )
    if is_inserted
      num_better += 1
      #print("candidate = "); show(candidate); println("")
      #print("index = "); show(index); println("")
      #print("target = "); show(de.population[index,:]); println("")
      old = de.population[:,index]
      de.population[:,index] = candidate
    end
  end
  num_better
end

# Now we can create specific DE optimizers that are commonly used in the
# literature.

# The most used DE/rand/1/bin.
function de_rand_1_bin(options = @compat Dict{Symbol,Any}(); sampler = random_sampler, name = "DE/rand/1/bin")
  opts = Parameters(options, DE_DefaultOptions)
  DiffEvoOpt(name, opts[:Population], opts,
        FixedDiffEvoParameters(opts, size(opts[:Population], 2)), sampler,
        NoMutation(), DiffEvoRandBin1(), RandomBound(opts[:SearchSpace]))
end

function de_rand_2_bin(options = @compat Dict{Symbol,Any}(); sampler = random_sampler, name = "DE/rand/2/bin")
  opts = Parameters(options, DE_DefaultOptions)
  DiffEvoOpt(name, opts[:Population], opts,
        FixedDiffEvoParameters(opts, size(opts[:Population], 2)), sampler,
        NoMutation(), DiffEvoRandBin2(), RandomBound(opts[:SearchSpace]))
end

# The most used DE/rand/1/bin with "local geography" via radius limited sampling.
function de_rand_1_bin_radiuslimited(options = @compat Dict{Symbol,Any}())
  de_rand_1_bin(options; sampler = radius_limited_sampler, name = "DE/rand/1/bin/radiuslimited")
end

function de_rand_2_bin_radiuslimited(options = @compat Dict{Symbol,Any}())
  de_rand_2_bin(options; sampler = radius_limited_sampler, name = "DE/rand/2/bin/radiuslimited")
end
