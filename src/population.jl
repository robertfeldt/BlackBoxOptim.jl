abstract Population
abstract PopulationWithFitness{F} <: Population

# The simplest population implementation -- a matrix of floats, each column is an individual.
typealias PopulationMatrix Matrix{Float64}

popsize(pop::PopulationMatrix) = size(pop, 2)
numdims(pop::PopulationMatrix) = size(pop, 1)
params_mean(pop::PopulationMatrix) = mean(pop, 1)
params_std(pop::PopulationMatrix) = std(pop, 1)

popsize{F}(pop::Vector{Candidate{F}}) = length(pop)
numdims{F}(pop::Vector{Candidate{F}}) = isempty(pop) ? 0 : length(pop[1].params)

type FitPopulation{F} <: PopulationWithFitness{F}
  # The population is a matrix of floats, each column being an individual.
  individuals::PopulationMatrix

  nafitness::F
  fitness::Vector{F}

  candi_pool::Vector{Candidate{F}} # pool of reusable candidates

  function FitPopulation(individuals::PopulationMatrix, nafitness::F, fitness::Vector{F})
    popsize(individuals) == length(fitness) || throw(DimensionMismatch("Fitness vector length does not match the population size"))
    new(individuals, nafitness, fitness, @compat(Vector{Candidate{F}}()))
  end
end

FitPopulation{F}(individuals::PopulationMatrix, nafitness::F,
                 fitnesses::Vector{F} = fill(nafitness, popsize(individuals))) =
  FitPopulation{F}(individuals, nafitness, fitnesses)

FitPopulation(fs::FitnessScheme, individuals::PopulationMatrix) =
  FitPopulation(individuals, nafitness(fs))

FitPopulation(fs::FitnessScheme, popsize::Int = 100, dims::Int = 1) =
  # FIXME use PopulationMatrix(), this is v0.3 workaround
  FitPopulation(fs, Array(Float64, dims, popsize))

popsize(pop::FitPopulation) = popsize(pop.individuals)
numdims(pop::FitPopulation) = numdims(pop.individuals)
params_mean(pop::FitPopulation) = mean(pop.individuals, 1)
params_std(pop::FitPopulation) = std(pop.individuals, 1)

fitness(pop::FitPopulation, ix::Int) = pop.fitness[ix]

Base.getindex(pop::FitPopulation, rows, cols) = pop.individuals[rows, cols]
Base.getindex(pop::FitPopulation, ::Colon, cols) = pop.individuals[:, cols] # FIXME remove v0.3 workaround
Base.getindex(pop::FitPopulation, indi_ixs) = pop.individuals[:, indi_ixs]

function Base.append!{F}(pop::FitPopulation{F}, extra_pop::FitPopulation{F})
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

# get unitialized individual from a pool, or create one, if it's empty
function acquire_candi{F}(pop::FitPopulation{F})
  if isempty(pop.candi_pool)
    return Candidate{F}(@compat(Vector{Float64}(numdims(pop))), -1, pop.nafitness)
  end
  res = pop!(pop.candi_pool)
  # reset reference to genetic operation
  res.op = NO_GEN_OP
  res.tag = 0
  return res
end

# get an individual from a pool and set it to ix-th individual from population
function acquire_candi(pop::FitPopulation, ix::Int)
    x = acquire_candi(pop)
    x.params[:] = pop[ix] # FIXME might be suboptimal until Julia has array refs
    x.index = ix
    x.fitness = fitness(pop, ix)
    return x
end

# get an individual from a pool and set it to another candidate
acquire_candi{F}(pop::FitPopulation{F}, candi::Candidate{F}) = copy!(acquire_candi(pop), candi)

# put the candidate back to the pool
release_candi{F}(pop::FitPopulation{F}, candi::Candidate{F}) = push!(pop.candi_pool, candi)

# put the candidate into the population
function accept_candi!{F}(pop::FitPopulation{F}, candi::Candidate{F})
  pop.individuals[:, candi.index] = candi.params
  pop.fitness[candi.index] = candi.fitness
  release_candi(pop, candi)
end

function reset_fitness!{F}(candi::Candidate{F}, pop::FitPopulation{F})
  candi.fitness = pop.nafitness
  return candi
end

candi_pool_size(pop::FitPopulation) = length(pop.candi_pool)

# default population generation
function population(problem::OptimizationProblem, options::Parameters = EMPTY_PARAMS)
  if !haskey(options, :Population)
      pop = rand_individuals_lhs(search_space(problem), get(options, :PopulationSize, 50))
  else
     pop = options[:Population]
  end
  if isa(pop, Population)
    return pop
  elseif isa(pop, PopulationMatrix)
    return FitPopulation(fitness_scheme(problem), pop)
  else
    throw(ArgumentError("\"Population\" parameter is of unsupported type: $(typeof(pop))"))
  end
end
