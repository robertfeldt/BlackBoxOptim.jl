abstract GeneticOperator
abstract MutationOperator <: GeneticOperator
abstract CrossoverOperator{NP,NC} <: GeneticOperator

apply{T <: Real}(o::MutationOperator, parents::Vector{Vector{T}}) = map(p -> apply(o, p), parents)

numchildren(o::GeneticOperator) = 1
numparents(o::MutationOperator) = 1 # But it will apply to each parent separately if given more than one...
numparents{NP,NC}(o::CrossoverOperator{NP,NC}) = NP
numchildren{NP,NC}(o::CrossoverOperator{NP,NC}) = NC

# mutation operator that does nothing
type NoMutation <: MutationOperator end
function apply!(mo::NoMutation, target) end

include("mutation/polynomial_mutation.jl")
include("mutation/mutation_clock.jl")
include("crossover/simulated_binary_crossover.jl")
include("crossover/differential_evolution_crossover.jl")
