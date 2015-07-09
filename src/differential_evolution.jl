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

type DiffEvoOpt{P<:DiffEvoParameters,S<:IndividualsSelector,M<:MutationOperator,X<:DiffEvoCrossoverOperator,E<:EmbeddingOperator} <: DifferentialEvolutionOpt
  # TODO when sampler and bound would be parameterized, name is no longer required -- as everything is seen in the type name
  name::ASCIIString

  # A population is a matrix of floats, individuals stored in columns.
  population::Array{Float64, 2}

  # Set of operators that together define a specific DE strategy.
  params::P        # adjust crossover parameters after fitness calculation
  select::S        # random individuals selector
  mutate::M        # mutation operator
  crossover::X     # crossover operator
  embed::E         # embedding operator
end

function DiffEvoOpt{P<:DiffEvoParameters, S<:IndividualsSelector,
                    M<:MutationOperator, X<:DiffEvoCrossoverOperator, E<:EmbeddingOperator}(
    name::ASCIIString, pop, params::P,
    select::S = S(), mutate::M = M(), crossover::X = X(), embed::E = E())
  DiffEvoOpt{P,S,M,X,E}(name, pop, params, select, mutate, crossover, embed)
end

popsize(pop::Matrix{Float64}) = size(pop,2)
popsize(opt::DifferentialEvolutionOpt) = popsize(population(opt))

make_candidate(pop::Matrix{Float64}, i::Int) = Candidate{Float64}(pop[:, i], i)

# Ask for a new candidate object to be evaluated, and a list of individuals
# it should be ranked with. The individuals are supplied as an array of tuples
# with the individual and its index.
function ask(de::DiffEvoOpt)
  # Sample parents and target
  indices = select(de.select, de.population, 1 + numparents(de.crossover))
  parent_indices = indices[1:numparents(de.crossover)]
  #println("parent_indices = $(parent_indices)")
  target_index = indices[end]
  target = make_candidate(de.population, target_index) # FIXME we can cache target somewhere
  #println("target[$(target_index)] = $(target)")

  # Crossover parents and target
  @assert numchildren(de.crossover)==1
  trial = copy(target) # FIXME reuse some trial vector
  apply!( de.crossover,
          crossover_parameters( de.params, target_index )...,
          trial.params, de.population, parent_indices )
  # embed the trial parameter vector into the search space
  apply!(de.embed, trial.params, de.population, [target_index])
  #println("trial = $(trial)")

  # Return the candidates that should be ranked as tuples including their
  # population indices.
  return Candidate{Float64}[trial, target]
end

# Tell the optimizer about the ranking of candidates. Returns the number of
# better candidates that were inserted into the population.
function tell!{F}(de::DiffEvoOpt,
  # archive::Archive, # Skip for now
  rankedCandidates::Vector{Candidate{F}})
  num_candidates = length(rankedCandidates)
  num_better = 0
  for i in 1:div(num_candidates, 2)
    candi = rankedCandidates[i]
    is_inserted = candi.params != de.population[:,candi.index]
    adjust!( de.params, candi.index, is_inserted )
    if is_inserted
      num_better += 1
      #print("candidate = "); show(candidate); println("")
      #print("index = "); show(index); println("")
      #print("target = "); show(de.population[index,:]); println("")
      #old = de.population[:,index]
      de.population[:,candi.index] = candi.params
    end
  end
  num_better
end

# Now we can create specific DE optimizers that are commonly used in the
# literature.

function diffevo(problem::OptimizationProblem, name::ASCIIString,
                 select::IndividualsSelector = SimpleSelector(),
                 crossover::DiffEvoCrossoverOperator = DiffEvoRandBin1(),
                 options = @compat Dict{Symbol,Any}())
  opts = Parameters(options, DE_DefaultOptions)
  pop = opts[:Population]
  DiffEvoOpt(name, pop, FixedDiffEvoParameters(opts, popsize(pop)), select,
        NoMutation(), crossover, RandomBound(search_space(problem)))
end

# The most used DE/rand/1/bin.
de_rand_1_bin(problem::OptimizationProblem,
              options = @compat(Dict{Symbol,Any}()),
              name = "DE/rand/1/bin") = diffevo(problem, name, SimpleSelector(), DiffEvoRandBin1(), options)

de_rand_2_bin(problem::OptimizationProblem,
              options = @compat(Dict{Symbol,Any}()),
              name = "DE/rand/2/bin") = diffevo(problem, name, SimpleSelector(), DiffEvoRandBin2(), options)

# The most used DE/rand/1/bin with "local geography" via radius limited sampling.
de_rand_1_bin_radiuslimited(problem::OptimizationProblem,
                            options = @compat(Dict{Symbol,Any}()),
                            name = "DE/rand/1/bin/radiuslimited") =
    diffevo(problem, name, RadiusLimitedSelector(Parameters(options, DE_DefaultOptions)[:SamplerRadius]),
            DiffEvoRandBin1(), options)

de_rand_2_bin_radiuslimited(problem::OptimizationProblem,
                            options = @compat(Dict{Symbol,Any}()),
                            name = "DE/rand/2/bin/radiuslimited") =
    diffevo(problem, name, RadiusLimitedSelector(Parameters(options, DE_DefaultOptions)[:SamplerRadius]),
            DiffEvoRandBin2(), options)
