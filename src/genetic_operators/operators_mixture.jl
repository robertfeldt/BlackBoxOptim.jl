# Randomly selects the mutator from the vector
# according to its weight and applies it
type FixedGeneticOperatorsMixture <: GeneticOperatorsMixture
    operators::Vector{GeneticOperator} # available operations
    weights::WeightVec{Float64}        # fixed weights

    function FixedGeneticOperatorsMixture(
        operators::Vector{GeneticOperator},
        rates::Vector{Float64} = fill(1.0/length(operators), length(operators)) # defaults to uniform distribution of rates
    )
      length(operators) == length(rates) || throw(DimensionMismatch("Number of mutators does not match the number of their rates"))
      new(operators, weights(rates))
    end
end

# default implementation of apply!() for operators mixture
function apply!{T<:Real}(opmix::GeneticOperatorsMixture, v::Vector{T})
  op, tag = next(opmix)
  apply!(op, v)
end

function next(opmix::FixedGeneticOperatorsMixture)
  i = sample(1:length(opmix.operators), opmix.weights)
  return opmix.operators[i], i
end
