# Randomly selects the mutator from the vector
# according to its weight and applies it
type MutationMixture <: MutationOperator
    mutators::Vector{MutationOperator} # available mutations
    weights::WeightVec{Float64}        # mutation weights

    function MutationMixture(mutators::Vector{MutationOperator}, rates::Vector{Float64})
      length(mutators) == length(rates) || throw(DimensionMismatch("Number of mutators does not match the number of their rates"))
      new(mutators, weights(rates))
    end
    MutationMixture(mutators::Vector{MutationOperator}) =
      # uniform distribution of rates of diff
      new(mutators, weights(fill(1.0/length(mutators), length(mutators))))
end

function apply!{T<:Real}(mm::MutationMixture, v::Vector{T})
  if isempty(mm.mutators) return v end # no change

  m = sample(mm.mutators, mm.weights) # randomly select the mutation to apply
  return apply!(m, v)
end
