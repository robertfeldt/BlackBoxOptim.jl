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
immutable NoMutation <: MutationOperator end
apply!(mo::NoMutation, target) = target

# placeholder for no-effect genetic operations
const NO_GEN_OP = NoMutation()

# adjust the internal parameters of the genetic operator
# default implementation does nothing
function adjust!{F}(op::GeneticOperator, tag::Int, indi_index::Int, new_fitness::F, old_fitness::F, is_improved::Bool) end

# a mixture of genetic operators,
# use next() to choose the next operator from the mixture
abstract GeneticOperatorsMixture <: GeneticOperator

include("operators_mixture.jl")
include("mutation/polynomial_mutation.jl")
include("mutation/mutation_clock.jl")
include("crossover/simulated_binary_crossover.jl")
include("crossover/differential_evolution_crossover.jl")
include("embedding/random_bound.jl")
include("selector/simple.jl")
include("selector/radius_limited.jl")
