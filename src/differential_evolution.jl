using Distributions

const DE_DefaultOptions = chain(DEX_DefaultOptions, ParamsDict(
    :SamplerRadius => 8,
))

# FIXME DifferentialEvolution is just a specific case of this optimizer,
# should it be called EvolutionaryOptimizer?
mutable struct DiffEvoOpt{P<:Population,S<:IndividualsSelector,M<:GeneticOperator,E<:EmbeddingOperator} <: PopulationOptimizer
    # TODO when sampler and bound would be parameterized, name is no longer required -- as everything is seen in the type name
    name::String

    population::P

    # Set of operators that together define a specific DE strategy.
    select::S        # random individuals selector
    modify::M        # genetic operator
    embed::E         # embedding operator

    function DiffEvoOpt(name::String, pop::P,
                        select::S = S(), modify::M = M(), embed::E = E()) where
            {P<:Population, S<:IndividualsSelector,
             M<:GeneticOperator, E<:EmbeddingOperator}
        new{P,S,M,E}(name, pop, select, modify, embed)
    end
end

# trace current optimization state,
# Called by OptRunController trace_progress()
function trace_state(io::IO, de::DiffEvoOpt, mode::Symbol)
    if mode != :compact
        println(io, "DE modify state:")
        trace_state(io, de.modify, mode)
    end
end

ask(de::DiffEvoOpt) = evolve(de, de.modify)

function evolve(de::DiffEvoOpt, op::GeneticOperatorsMixture)
    sel_op, tag = next(op)
    candidates = evolve(de, sel_op) # use the selected operator of the mixture
    # override what was set by sel_op so that adjust!() gets called for the mixture op
    for candi in candidates
        candi.extra = op
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

"""
Post-process the evolved pair.
"""
function evolved_pair(de::DiffEvoOpt, target::Candidate{F}, trial::Candidate{F},
                      op::GeneticOperator, tag::Int = 0) where F
    # embed the trial parameter vector into the search space
    apply!(de.embed, trial.params, de.population, target.index)
    target.extra = trial.extra = op
    target.tag = trial.tag = tag
    if trial.params != target.params
        reset_fitness!(trial, de.population)
    end

    return Candidate{F}[trial, target]
end

function tell!(de::DiffEvoOpt,
        # archive::Archive, # Skip for now
        rankedCandidates::Vector{Candidate{F}}) where F
    n_acceptable_candidates = length(rankedCandidates)÷2
    num_better = 0
    for i in eachindex(rankedCandidates)
        candi = rankedCandidates[i]
        # accept the modified individuals from the top ranked half
        if i <= n_acceptable_candidates
            is_improved = candi.params != viewer(population(de), candi.index)
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

"""
Adjust the parameters of the method after the candidate evaluation.
"""
function adjust!(de::DiffEvoOpt, candi::Candidate, is_improved::Bool)
    # adjust the parameters of the operation
    old_fitness = fitness(population(de), candi.index)
    if isnan(old_fitness)
        old_fitness = candi.fitness
    end

    adjust!(candi.extra::GeneticOperator, candi.tag, candi.index, candi.fitness,
            fitness(population(de), candi.index), is_improved)
end

# Now we can create specific DE optimizers that are commonly used in the
# literature.

function diffevo(problem::OptimizationProblem, options::Parameters, name::String,
                 select::IndividualsSelector = SimpleSelector(),
                 crossover::DiffEvoCrossoverOperator = DiffEvoRandBin1(chain(DE_DefaultOptions, options)))
    opts = chain(DE_DefaultOptions, options)
    pop = population(problem, opts)
    DiffEvoOpt(name, pop, select, crossover,
               RandomBound(search_space(problem)))
end

"""
The most used `DE/rand/1/bin` variant of differential evolution.
"""
de_rand_1_bin(problem::OptimizationProblem,
              options::Parameters = EMPTY_PARAMS,
              name = "DE/rand/1/bin") = diffevo(problem, options, name)

de_rand_2_bin(problem::OptimizationProblem,
              options::Parameters = EMPTY_PARAMS,
              name = "DE/rand/2/bin") = diffevo(problem, options, name,
                                                SimpleSelector(),
                                                DiffEvoRandBin2(chain(DE_DefaultOptions, options)))

"""
The most used `DE/rand/1/bin` variant with "local geography" via radius-limited sampling.
"""
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
            DiffEvoRandBin2(chain(DE_DefaultOptions, options)))
