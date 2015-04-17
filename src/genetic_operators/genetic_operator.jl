abstract GeneticOperator
abstract MutationOperator <: GeneticOperator
abstract CrossoverOperator <: GeneticOperator

# Apply a single-value crossover on all values of parent vectors (p1 & p2) to get
# two child vectors.
function apply{T <: Real}(xo::CrossoverOperator, p1::Vector{T}, p2::Vector{T}, c1::Vector{T}, c2::Vector{T})
  for i in 1:length(p1)
    n1, n2 = apply(xo, p1[i], p2[i])
    c1[i] = n1
    c2[i] = n2
  end
  return c1, c2
end

apply!{T <: Real}(xo::CrossoverOperator, p1::Vector{T}, p2::Vector{T}) = apply(xo, p1, p2, p1, p2)

include("mutation/polynomial_mutation.jl")
include("mutation/mutation_clock.jl")
include("crossover/simulated_binary_crossover.jl")
