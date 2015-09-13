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
function apply!{T<:Real}(opmix::GeneticOperatorsMixture, v::Vector{T}, target_index::Int)
  op, tag = next(opmix)
  apply!(op, v, target_index)
end

function next(opmix::FixedGeneticOperatorsMixture)
  i = sample(1:length(opmix.operators), opmix.weights)
  return opmix.operators[i], i
end

# Frequency-adapting genetic operators mixture
type FAGeneticOperatorsMixture <: GeneticOperatorsMixture
    operators::Vector{GeneticOperator} # available operators
    fa::FrequencyAdapter               # adapter of operator frequencies

    function FAGeneticOperatorsMixture(operators::Vector{GeneticOperator}; c = 1E-2, eta = 1E-2, pmin = 0.01, pmax = 1.0)
      new(operators, FrequencyAdapter(length(operators); c = c, eta = eta, pmin = pmin, pmax = pmax))
    end
end

frequencies(opmix::FAGeneticOperatorsMixture) = frequencies(opmix.fa)

function next(opmix::FAGeneticOperatorsMixture)
  i = next(opmix.fa)
  return opmix.operators[i], i
end

function adjust!{F}(opmix::FAGeneticOperatorsMixture, op_index::Int, candi_index::Int, new_fitness::F, old_fitness::F, is_improved::Bool)
  # KLUDGE we don't know the fitness scheme, but if there is no improvement in fitness,
  # new_fitness==old_fitness, so it's ok to take the abs()
  update!(opmix.fa, op_index, is_improved ? abs(new_fitness-old_fitness) : 0.0)

  # also adjust the actual operator
  adjust!(opmix.operators[op_index], 0, candi_index, new_fitness, old_fitness, is_improved)
end

function trace_state(io::IO, fa::FAGeneticOperatorsMixture)
  println(io, "FrequencyAdapting mixture rates:")
  for i in eachindex(fa.operators)
    @printf(io, " %s=%.2f", typeof(fa.operators[i]), fa.fa.p[i])
  end
  println(io)
end
