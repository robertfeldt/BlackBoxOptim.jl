"""
  Abstract genetic operator that transforms individuals in the population.
"""
abstract GeneticOperator

"""
  Modifies (mutates) one individual.

  The concrete implementations must provide `apply!()` method.
"""
abstract MutationOperator <: GeneticOperator

"""
  Modifies `NC` "children" by transferring some information from `NP` "parents".

  The concrete implementations must provide `apply!()` method.
"""
abstract CrossoverOperator{NP,NC} <: GeneticOperator

"""
  Embeds(projects) the individual into the search space.

  The concrete implementations must provide `apply!()` method.
"""
abstract EmbeddingOperator <: GeneticOperator

"""
  Selects the individuals from the population.

  The concrete implementations must provide `select()` method.
"""
abstract IndividualsSelector

"""
  `select(selector<:IndividualsSelector, population, numSamples::Int)`

  Select `numSamples` random candidates from the `population`.
"""
function select(::IndividualsSelector, population, numSamples::Int) end

apply{T <: Real}(o::MutationOperator, parents::Vector{Vector{T}}) = map(p -> apply(o, p), parents)

numchildren(o::GeneticOperator) = 1
numparents(o::MutationOperator) = 1 # But it will apply to each parent separately if given more than one...

numparents{NP,NC}(o::CrossoverOperator{NP,NC}) = NP
numchildren{NP,NC}(o::CrossoverOperator{NP,NC}) = NC

numparents(o::EmbeddingOperator) = 1
numchildren(o::EmbeddingOperator) = 1

# wrapper for multi-children variant of apply!() for single-child xover operators
function apply!{NP,T<:AbstractVector{Float64}}(xover::CrossoverOperator{NP,1}, targets::AbstractVector{T}, target_indices::AbstractVector{Int}, pop, parentIndices)
    length(targets) == length(target_indices) || throw(ArgumentError("The number of target doesn't match the number of their indices"))
    for i in eachindex(target_indices)
        apply!(xover, targets[i], target_indices[i], pop, parentIndices)
    end
    targets
end

"""
  `MutationOperator` that does nothing.
"""
immutable NoMutation <: MutationOperator end
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
function adjust!{F}(op::GeneticOperator, tag::Int, indi_index::Int, new_fitness::F, old_fitness::F, is_improved::Bool) end

"""
  `trace_state(io, op::GeneticOperator, mode::Symbol)`

  Trace the state of the operator.
  Called by `trace_progress()` during `OptRunController` run by some of the genetic optimizers.

  Override the method to trace the state of your genetic operator.
"""
function trace_state(io::IO, op::GeneticOperator, mode::Symbol) end

"""
  A mixture of genetic operators,
  use `next()` to choose the next operator from the mixture.
"""
abstract GeneticOperatorsMixture <: GeneticOperator

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
