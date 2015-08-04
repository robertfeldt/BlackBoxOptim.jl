using Distributions

const DE_DefaultOptions = @compat Dict{Symbol,Any}(
  :f => 0.6,
  :cr => 0.7,
  :SamplerRadius => 8,
)

# FIXME DifferentialEvolution is just a specific case of this optimizer,
# should it be called EvolutionaryOptimizer?
type DiffEvoOpt{P<:Population,S<:IndividualsSelector,M<:GeneticOperator,E<:EmbeddingOperator} <: PopulationOptimizer
  # TODO when sampler and bound would be parameterized, name is no longer required -- as everything is seen in the type name
  name::ASCIIString

  population::P

  # Set of operators that together define a specific DE strategy.
  select::S        # random individuals selector
  modify::M        # genetic operator
  embed::E         # embedding operator
end

function DiffEvoOpt{P<:Population, S<:IndividualsSelector,
                    M<:GeneticOperator, E<:EmbeddingOperator}(
    name::ASCIIString, pop::P,
    select::S = S(), modify::M = M(), embed::E = E())
  DiffEvoOpt{P,S,M,E}(name, pop, select, modify, embed)
end

# Ask for a new candidate object to be evaluated, and a list of individuals
# it should be ranked with. The individuals are supplied as an array of tuples
# with the individual and its index.
function ask(de::DiffEvoOpt)
  evolve(de, de.modify)
end

function evolve(de::DiffEvoOpt, op::GeneticOperatorsMixture)
  sel_op, tag = next(op)
  candidates = evolve(de, sel_op) # use the selected operator of the mixture
  # override what was set by sel_op so that adjust!() gets called for the mixture op
  for candi in candidates
    candi.op = op
    candi.tag = tag
  end
  return candidates
end

function evolve(de::DiffEvoOpt, mutate::MutationOperator)
  target_index = select(de.select, de.population, 1)[1]
  target = acquire_candi(de.population, target_index)
  trial = acquire_candi(de.population, target)
  apply!(mutate, trial.params, target_index)
  return evolved_pair(de, target, trial, mutate, 0)
end

function evolve(de::DiffEvoOpt, crossover::CrossoverOperator)
  # Sample parents and target
  indices = select(de.select, de.population, 1 + numparents(crossover))
  parent_indices = indices[1:numparents(crossover)]
  #println("parent_indices = $(parent_indices)")
  target_index = indices[end]
  target = acquire_candi(de.population, target_index)
  #println("target[$(target_index)] = $(target)")

  # Crossover parents and target
  @assert numchildren(crossover)==1
  trial = acquire_candi(de.population, target)
  apply!(crossover, trial.params, target_index,
         de.population, parent_indices)
  return evolved_pair(de, target, trial, crossover, 0)
end

# post processing of evolved pair
function evolved_pair{F}(de::DiffEvoOpt, target::Candidate{F}, trial::Candidate{F}, op::GeneticOperator, tag::Int = 0)
  # embed the trial parameter vector into the search space
  apply!(de.embed, trial.params, de.population, [target.index])
  target.op = trial.op = op
  target.tag = trial.tag = 0
  if trial.params != target.params
    reset_fitness!(trial, de.population)
  end

  return Candidate{F}[trial, target]
end

# Tell the optimizer about the ranking of candidates. Returns the number of
# better candidates that were inserted into the population.
function tell!{F}(de::DiffEvoOpt,
  # archive::Archive, # Skip for now
  rankedCandidates::Vector{Candidate{F}})
  n_acceptable_candidates = length(rankedCandidates)÷2
  num_better = 0
  for i in eachindex(rankedCandidates)
    candi = rankedCandidates[i]
    # accept the modified individuals from the top ranked half
    if i <= n_acceptable_candidates
      is_improved = candi.params != population(de)[candi.index]
      adjust!(de, candi, is_improved)
    else
      is_improved = false
    end
    if is_improved
      num_better += 1
      #print("candidate = "); show(candidate); println("")
      #print("index = "); show(index); println("")
      #print("target = "); show(de.population[index,:]); println("")
      #old = de.population[:,index]
      accept_candi!(de.population, candi)
    else
      # just return candidate to the pool
      release_candi(de.population, candi)
    end
  end
  num_better
end

# adjust the parameters of the method
function adjust!(de::DiffEvoOpt, candi::Candidate, is_improved::Bool)
  # adjust the parameters of the operation
  old_fitness = fitness(population(de), candi.index)
  if isnan(old_fitness)
    old_fitness = candi.fitness
  end

  adjust!(candi.op, candi.tag, candi.index, candi.fitness,
          fitness(population(de), candi.index), is_improved)
end

# Now we can create specific DE optimizers that are commonly used in the
# literature.

function diffevo(problem::OptimizationProblem, options::Parameters, name::ASCIIString,
                 select::IndividualsSelector = SimpleSelector(),
                 crossover::DiffEvoCrossoverOperator = convert(DiffEvoRandBin1, chain(DE_DefaultOptions, options)))
  opts = chain(DE_DefaultOptions, options)
  pop = population(problem, opts)
  DiffEvoOpt(name, pop, select, crossover,
             RandomBound(search_space(problem)))
end

# The most used DE/rand/1/bin.
de_rand_1_bin(problem::OptimizationProblem,
              options::Parameters = EMPTY_PARAMS,
              name = "DE/rand/1/bin") = diffevo(problem, options, name)

de_rand_2_bin(problem::OptimizationProblem,
              options::Parameters = EMPTY_PARAMS,
              name = "DE/rand/2/bin") = diffevo(problem, options, name,
                                                SimpleSelector(),
                                                convert(DiffEvoRandBin2, chain(DE_DefaultOptions, options)))

# The most used DE/rand/1/bin with "local geography" via radius limited sampling.
de_rand_1_bin_radiuslimited(problem::OptimizationProblem,
                            options::Parameters = EMPTY_PARAMS,
                            name = "DE/rand/1/bin/radiuslimited") =
    diffevo(problem, options, name,
            RadiusLimitedSelector(chain(DE_DefaultOptions, options)[:SamplerRadius]))

de_rand_2_bin_radiuslimited(problem::OptimizationProblem,
                            options::Parameters = EMPTY_PARAMS,
                            name = "DE/rand/2/bin/radiuslimited") =
    diffevo(problem, options, name,
            RadiusLimitedSelector(chain(DE_DefaultOptions, options)[:SamplerRadius]),
            convert(DiffEvoRandBin2, chain(DE_DefaultOptions, options)))
