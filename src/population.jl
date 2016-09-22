"""
  The base abstract type for the collection of candidate solutions
  in the population-based optimization methods.
"""
abstract Population
"""
  The base abstract types for population that also stores the candidates
  fitness.

  `F` is the fitness type.
"""
abstract PopulationWithFitness{F} <: Population

"""
  The simplest `Population` implementation -- a matrix of floats, each column is an individual.
"""
typealias PopulationMatrix Matrix{Float64}

popsize(pop::PopulationMatrix) = size(pop, 2)
numdims(pop::PopulationMatrix) = size(pop, 1)
params_mean(pop::PopulationMatrix) = mean(pop, 1)
params_std(pop::PopulationMatrix) = std(pop, 1)

popsize{F}(pop::Vector{Candidate{F}}) = length(pop)
numdims{F}(pop::Vector{Candidate{F}}) = isempty(pop) ? 0 : length(pop[1].params)

viewer(pop::PopulationMatrix, indi_ix) = view(pop, :, indi_ix)

"""
  The default implementation of `PopulationWithFitness{F}`.
"""
type FitPopulation{F} <: PopulationWithFitness{F}
  # The population is a matrix of floats, each column being an individual.
  individuals::PopulationMatrix

  nafitness::F
  fitness::Vector{F}
  ntransient::Int                  # how many transient members are in the population

  candi_pool::Vector{Candidate{F}} # pool of reusable candidates

  function FitPopulation(individuals::PopulationMatrix, nafitness::F, fitness::Vector{F}; ntransient::Integer=0)
    popsize(individuals) == length(fitness) || throw(DimensionMismatch("Fitness vector length does not match the population size"))
    new(individuals, nafitness, fitness, ntransient, Vector{Candidate{F}}())
  end
end

FitPopulation{F}(individuals::PopulationMatrix, nafitness::F,
                 fitnesses::Vector{F} = fill(nafitness, popsize(individuals));
                 ntransient::Integer=0) =
  FitPopulation{F}(individuals, nafitness, fitnesses, ntransient=ntransient)

FitPopulation(fs::FitnessScheme, individuals::PopulationMatrix;
              ntransient::Integer=0) =
  FitPopulation(individuals, nafitness(fs), ntransient=ntransient)

FitPopulation(fs::FitnessScheme, popsize::Int = 100, dims::Int = 1;
              ntransient::Integer=0) =
  FitPopulation(fs, PopulationMatrix(dims, popsize); ntransient=ntransient)

popsize(pop::FitPopulation) = popsize(pop.individuals)-pop.ntransient
numdims(pop::FitPopulation) = numdims(pop.individuals)

# resize the population
function Base.resize!(pop::FitPopulation, newpopsize::Integer)
    new_individuals = PopulationMatrix(numdims(pop), newpopsize+pop.ntransient)
    new_individuals[:, 1:min(newpopsize,popsize(pop))] = view(pop.individuals, :, 1:min(newpopsize,popsize(pop)))
    pop.individuals = new_individuals
    resize!(pop.fitness, newpopsize + pop.ntransient)
    pop
end

# indices of the persistent individuals
@inline persistent_range(pop::FitPopulation) = 1:popsize(pop)
# indices of the transient individuals
@inline transient_range(pop::FitPopulation) = popsize(pop):(popsize(pop)+pop.ntransient-1)
persistent_individuals(pop::FitPopulation) = view(pop.individuals, :, persistent_range(pop))

params_mean(pop::FitPopulation) = mean(persistent_individuals(pop), 1)
params_std(pop::FitPopulation) = std(persistent_individuals(pop), 1)

fitness(pop::FitPopulation, ix::Int) = pop.fitness[ix]

Base.getindex(pop::FitPopulation, rows, cols) = pop.individuals[rows, cols]
Base.getindex(pop::FitPopulation, indi_ixs) = pop.individuals[:, indi_ixs]

"""
    viewer(population, individual_index)

    Get vector-slice of the population matrix for the specified
    individual.
    Does not allocate any additional space, while still providing the same lookup performance.
"""
viewer(pop::FitPopulation, indi_ix) = view(pop.individuals, :, indi_ix)

Base.setindex!(pop::FitPopulation, indi::Individual, ::Colon, indi_ix::Integer) =
    setindex!(pop, indi, indi_ix)

function Base.setindex!(pop::FitPopulation, indi::Individual, indi_ix::Integer)
    pop.individuals[:, indi_ix] = indi
    pop.fitness[indi_ix] = pop.nafitness
    pop
end

function Base.setindex!{F}(pop::FitPopulation{F}, indi::FitIndividual{F}, indi_ix::Integer)
    pop.individuals[:, indi_ix] = params(indi)
    pop.fitness[indi_ix] = fitness(indi)
    pop
end

function Base.append!{F}(pop::FitPopulation{F}, extra_pop::FitPopulation{F})
  pop.ntransient == 0 || throw(error("Appending to the population with transients not supported (yet)"))
  numdims(pop) == numdims(extra_pop) ||
    throw(DimensionMismatch("Cannot append population, "*
                            "the number of parameters differs "*
                            "($(numdims(pop)) vs $(numdims(extra_pop)))"))
  pop.individuals = hcat(pop.individuals, extra_pop.individuals)
  append!(pop.fitness, extra_pop.fitness)
  return pop
end

fitness_type{F}(pop::FitPopulation{F}) = F
candidate_type{F}(pop::FitPopulation{F}) = Candidate{F}

"""
  `acquire_candi(pop::FitPopulation[, {ix::Int, candi::Candidate}])`

  Get individual from a pool, or create one if the pool is empty.
  By default the individual is not initialized, but if `ix` or `candi` is specified,
  the corresponding fields of the new candidate are set to the given values.
"""
function acquire_candi{F}(pop::FitPopulation{F})
  if isempty(pop.candi_pool)
    return Candidate{F}(Individual(numdims(pop)), -1, pop.nafitness)
  end
  res = pop!(pop.candi_pool)
  # reset reference to genetic operation
  res.extra = NO_GEN_OP
  res.tag = 0
  return res
end

acquire_candis{F}(pop::FitPopulation{F}, n::Integer) =
    Candidate{F}[acquire_candi(pop) for _ in 1:n]

# Get an individual from a pool and sets it to ix-th individual from population.
function acquire_candi(pop::FitPopulation, ix::Int)
    x = acquire_candi(pop)
    copy!(x.params, viewer(pop, ix))
    x.index = ix
    x.fitness = fitness(pop, ix)
    return x
end

# Get an individual from a pool and sets it to another candidate.
acquire_candi{F}(pop::FitPopulation{F}, candi::Candidate{F}) = copy!(acquire_candi(pop), candi)

"""
  Put the candidate back to the pool.
"""
release_candi{F}(pop::FitPopulation{F}, candi::Candidate{F}) = push!(pop.candi_pool, candi)

"""
  Put the candidate back into the pool and copy the values
  into the corresponding individual of the population (`candi.index` should be set).
"""
function accept_candi!{F}(pop::FitPopulation{F}, candi::Candidate{F})
  pop.individuals[:, candi.index] = candi.params
  pop.fitness[candi.index] = candi.fitness
  release_candi(pop, candi)
end

"""
  Reset the candidate fitness.

  Need it when the candidate parameters have changed, but the stored fitness
  is still for the old parameter set.
"""
function reset_fitness!{F}(candi::Candidate{F}, pop::FitPopulation{F})
  candi.fitness = pop.nafitness
  return candi
end

candi_pool_size(pop::FitPopulation) = length(pop.candi_pool)

"""
  Generate a population for a given problem.

  The default method to generate a population, uses Latin Hypercube Sampling.
"""
function population{F}(problem::OptimizationProblem,
                       options::Parameters = EMPTY_PARAMS,
                       nafitness::F = nafitness(fitness_scheme(problem));
                       ntransient::Integer = 0)
  if !haskey(options, :Population)
      pop = rand_individuals_lhs(search_space(problem), get(options, :PopulationSize, 50) + ntransient)
  else
     pop = options[:Population]
  end
  if isa(pop, Population)
    return pop
  elseif isa(pop, PopulationMatrix)
    return FitPopulation(pop, nafitness, ntransient=ntransient)
  else
    throw(ArgumentError("\"Population\" parameter is of unsupported type: $(typeof(pop))"))
  end
end
