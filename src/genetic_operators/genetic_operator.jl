abstract GeneticOperator
abstract MutationOperator <: GeneticOperator
abstract CrossoverOperator <: GeneticOperator

apply{T <: Real}(o::MutationOperator, parents::Vector{Vector{T}}) = map(p -> apply(o, p), parents)

numchildren(o::GeneticOperator) = 1
numparents(o::MutationOperator) = 1 # But it will apply to each parent separately if given more than one...

abstract XoverOp2to2 <: CrossoverOperator
numparents(o::XoverOp2to2) = 2
numchildren(o::XoverOp2to2) = 2

abstract XoverOp3to1 <: CrossoverOperator
numparents(o::XoverOp3to1) = 3

abstract XoverOp5to1 <: CrossoverOperator
numparents(o::XoverOp5to1) = 5
numchildren(o::XoverOp5to1) = 1

apply{T <: Real}(xo::CrossoverOperator, parents::Vector{Vector{T}}) = apply(xo, parents...)

include("mutation/polynomial_mutation.jl")
include("mutation/mutation_clock.jl")
include("crossover/simulated_binary_crossover.jl")
include("crossover/differential_evolution_crossover.jl")
