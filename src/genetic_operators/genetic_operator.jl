"""
Abstract genetic operator that transforms individuals in the population.
"""
abstract type GeneticOperator end

"""
Modifies (mutates) one individual.

The concrete implementations must provide `apply!()` method.
"""
abstract type MutationOperator <: GeneticOperator end

"""
Modifies `NC` "children" by transferring some information from `NP` "parents".

The concrete implementations must provide `apply!()` method.
"""
abstract type CrossoverOperator{NP,NC} <: GeneticOperator end

"""
Embeds(projects) the individual into the search space.

The concrete implementations must provide `apply!()` method.
"""
abstract type EmbeddingOperator <: GeneticOperator end

"""
Selects the individuals from the population.

The concrete implementations must provide `select()` method.
"""
abstract type IndividualsSelector end

"""
    select(selector<:IndividualsSelector, population, numSamples::Int)

Select `numSamples` random candidates from the `population`.
"""
function select(::IndividualsSelector, population, numSamples::Int) end

apply(o::MutationOperator, parents::AbstractVector{<:AbstractVector{<:Real}}) =
    map(p -> apply(o, p), parents)

numchildren(o::GeneticOperator) = 1
numparents(o::MutationOperator) = 1 # But it will apply to each parent separately if given more than one...

numparents(o::CrossoverOperator{NP,NC}) where {NP,NC} = NP::Int
numchildren(o::CrossoverOperator{NP,NC}) where {NP,NC} = NC::Int

numparents(o::EmbeddingOperator) = 1
numchildren(o::EmbeddingOperator) = 1

# wrapper for multi-children variant of apply!() for single-child xover operators
function apply!(xover::CrossoverOperator{NP, 1},
                targets::AbstractVector{<:AbstractIndividual}, target_indices::AbstractVector{Int},
                pop, parentIndices) where NP
    length(targets) == length(target_indices) ||
        throw(ArgumentError("The number of target doesn't match the number of their indices"))
    for i in eachindex(target_indices)
        apply!(xover, targets[i], target_indices[i], pop, parentIndices)
    end
    targets
end

"""
`MutationOperator` that does nothing.
"""
struct NoMutation <: MutationOperator end
apply!(mo::NoMutation, target, target_index) = target

"""
Placeholder for no-effect genetic operations.
"""
const NO_GEN_OP = NoMutation()

"""
Adjust the internal parameters of the genetic operator `op` taking into account
the fitness change.

The default implementation does nothing.
"""
function adjust!(op::GeneticOperator, tag::Int, indi_index::Int,
                 new_fitness::F, old_fitness::F, is_improved::Bool) where F end

"""
    trace_state(io, op::GeneticOperator, mode::Symbol)

Trace the state of the operator.
Called by `trace_progress()` during `OptRunController` run by some of the genetic optimizers.

Override the method to trace the state of your genetic operator.
"""
function trace_state(io::IO, op::GeneticOperator, mode::Symbol) end

"""
A mixture of genetic operators,
use `next()` to choose the next operator from the mixture.
"""
abstract type GeneticOperatorsMixture <: GeneticOperator end

include("operators_mixture.jl")

include("mutation/mutation_clock.jl")
include("mutation/polynomial_mutation.jl")

include("crossover/mutation_wrapper.jl")
include("crossover/simulated_binary_crossover.jl")
include("crossover/simplex_crossover.jl")
include("crossover/differential_evolution_crossover.jl")
include("crossover/parent_centric_crossover.jl")
include("crossover/unimodal_normal_distribution_crossover.jl")

include("embedding/random_bound.jl")

include("selector/simple.jl")
include("selector/radius_limited.jl")
include("selector/tournament.jl")
