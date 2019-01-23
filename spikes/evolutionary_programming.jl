# This implements the common parts of Evolutionary Programming, i.e. the 
# population-based "culling" between generations. Different variants of EP
# can then be implemented by changing the mutation-operator etc.

PopulationBasedEvolutionaryProgrammingDefaultParameters = {
  :PopulationSize => 50,
  :InitialStandardDeviation => 3.0,  # Initial standard deviation
} 

# The specific functions of a particular EP strategy is attached to a strategy
# datum while the generic population-based EP is kept in a single type. This
# way we can easily change the strategy without having to rewrite the common parts
# throughout.
abstract type EvolutionaryProgrammingStrategy end

mutable struct PopulationBasedEvolutionaryProgramming <: PopulationOptimizer
  parameters::Parameters

  generation::Int
  dimension::Int
  population_size::Int

  # Pre-calc the tau values so we do not have to recalc it again and again during evolution.
  tau::Float64
  tauprim::Float64

  # Two populations that we switch between as being the current population.
  # The generation mod(number, 2) determines which one is current. This way
  # we can avoid continues memory allocation and freeing of the population 
  # structures.
  pop0::Array{Float64, 2}
  pop1::Array{Float64, 2}

  # Child population used in creating the next population
  children::Array{Float64, 2}

  # And the sigma values for each individual are kept in their separate matrices.
  # Exact same structure as for the populations and children matrices.
  sigmas0::Array{Float64, 2}
  sigmas1::Array{Float64, 2}
  childsigmas::Array{Float64, 2}

  PopulationBasedEvolutionaryProgramming(parameters) = begin

    dim = parameters[:NumDimensions]
    popsize = parameters[:PopulationSize]
    initial_population = rand_individuals_lhs(parameters[:SearchSpace], popsize)
    initial_sigmas = parameters[:InitialStandardDeviation] * ones(dims, popsize)

    new(parameters, 0, dims, popsize,
      1.0 / sqrt(2 * sqrt(popsize)), 1.0 / sqrt(2 * popsize),
      initial_population, zeros(Float64, dims, popsize), zeros(Float64, dims, popsize),
      initial_sigmas, zeros(Float64, dims, popsize), zeros(Float64, dims, popsize))

  end
end

function get_current_population(pep::PopulationBasedEvolutionaryProgramming)
  if mod(pep.generation, 2) == 0
    return pep.pop0, pep.sigmas0
  else
    return pep.pop1, pep.sigmas1
  end
end

# The generic EP strategy functions have the abstract type and thus specific
# strategies can fall back on them simply by not "overriding" that function.
# This implements a classic Gaussian EP with gaussian mutations (GEP).
has_ask_tell_interface(ep::EvolutionaryProgrammingStrategy) = false

function step(ep::EvolutionaryProgrammingStrategy)

  pep = ep.pep
  children = pep.children
  childsigmas = pep.childsigmas
  currentpop, currentsigmas, nextpop, nextsigmas = get_current_and_next_population_and_sigma(pep)

  # Mutate each parent to create a child.
  for i in 1:pep.population_size
    children[:,i], childsigmas[:,i] = mutate_individual(ep, currentpop, currentsigmas, 
      i, pep.dimension, pep.tauprim, pep.tau)
  end

  # Select the individuals to survive and copy them to the next population.
  select_winning_individuals(ep, currentpop, currentsigmas, children, childsigmas, nextpop, nextsigmas)

end

function select_winning_individuals(ep::EvolutionaryProgrammingStrategy,
  currentpop::Array{Float64, 2}, currentsigmas::Array{Float64, 2},
  children::Array{Float64, 2}, childsigmas::Array{Float64, 2},
  nextpop::Array{Float64, 2}, nextsigmas::Array{Float64, 2})

  

end

function mutate_individual(ep::EvolutionaryProgrammingStrategy, 
  population::Array{Float64, 2}, sigmas::Array{Float64, 2}, 
  i::Int, dims::Int, tauprim::Float64, tau::Float64)

  tauprimconst = tauprim * sample_tauprimconst_distribution(ep)

  childsigmas = sigmas[:,i] .* exp( tauprimconst + tau * sample_sigma_distribution(ep, dims)' )
  child = population[:,i] .+ childsigmas .* sample_mutation_distribution(ep, dims)'

  return child, childsigmas

end

function sample_tauprimconst_distribution(ep::EvolutionaryProgrammingStrategy)
  randn() # Standard is to use a N(0, 1) gaussian...
end

# The standard distribution is a Gaussian(0, 1)
function sample_mutation_distribution(ep::EvolutionaryProgrammingStrategy, numsamples::Int)
  randn(numsamples) # Standard is to use a N(0, 1) gaussians...
end

# The standard distribution for mutating the sigmas is a Gaussian(0, 1) for each dimension.
function sample_sigma_distribution(ep::EvolutionaryProgrammingStrategy, numsamples::Int)
  randn(numsamples) # Standard is to use a N(0, 1) gaussians...
end

# Specific EP variants are implemented in their own types and have the PBEP
# structure in them

mutable struct GaussianEvolutionaryProgramming <: EvolutionaryProgrammingStrategy
  pep::PopulationBasedEvolutionaryProgramming
end

