# abstract genetic operator that transforms individuals in the population
abstract GeneticOperator
# modifies one individual
abstract MutationOperator <: GeneticOperator
# modifies NC "children" by transferring some information from NP "parents"
abstract CrossoverOperator{NP,NC} <: GeneticOperator
# embeds(projects) the individual into the search space
abstract EmbeddingOperator <: GeneticOperator

# selects the individuals from the population
abstract IndividualsSelector

apply{T <: Real}(o::MutationOperator, parents::Vector{Vector{T}}) = map(p -> apply(o, p), parents)

numchildren(o::GeneticOperator) = 1
numparents(o::MutationOperator) = 1 # But it will apply to each parent separately if given more than one...

numparents{NP,NC}(o::CrossoverOperator{NP,NC}) = NP
numchildren{NP,NC}(o::CrossoverOperator{NP,NC}) = NC

numparents(o::EmbeddingOperator) = 1
numchildren(o::EmbeddingOperator) = 1

# mutation operator that does nothing
type NoMutation <: MutationOperator end
function apply!(mo::NoMutation, target) end

include("mutation/polynomial_mutation.jl")
include("mutation/mutation_clock.jl")
include("mutation/mutation_mixture.jl")
include("crossover/simulated_binary_crossover.jl")
include("crossover/differential_evolution_crossover.jl")
include("embedding/random_bound.jl")
include("selector/simple.jl")
include("selector/radius_limited.jl")
