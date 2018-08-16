"""
Randomly selects the genetic operator from the vector
according to its weight and applies it.
"""
struct FixedGeneticOperatorsMixture <: GeneticOperatorsMixture
    operators::Vector{GeneticOperator} # available operations
    weights::Weights{Float64, Float64, Vector{Float64}}  # fixed weights

    function FixedGeneticOperatorsMixture(
        operators::AbstractVector{<:GeneticOperator},
        rates::AbstractVector{Float64} = fill(1.0/length(operators), length(operators)) # defaults to uniform distribution of rates
    )
        length(operators) == length(rates) ||
            throw(DimensionMismatch("Number of mutators does not match the number of their rates"))
        new(GeneticOperator[op for op in operators], weights(rates))
    end
end

"""
Default implementation of `apply!()` for operators mixture.
"""
function apply!(opmix::GeneticOperatorsMixture,
                v::AbstractVector{<:Real}, target_index::Int)
    op, tag = next(opmix)
    apply!(op, v, target_index)
end

"""
    next(fa::FixedGeneticOperatorsMixture)

Gets the random genetic operator from the mixture.

The probability to select an operator is proportional to its weight.

Returns a tuple of the genetic operator and its index in the mix.
"""
function next(opmix::FixedGeneticOperatorsMixture)
    i = sample(1:length(opmix.operators), opmix.weights)
    return opmix.operators[i], i
end

"""
Frequency-adapting genetic operators mixture.
"""
struct FAGeneticOperatorsMixture <: GeneticOperatorsMixture
    operators::Vector{GeneticOperator} # available operators
    fa::FrequencyAdapter               # adapter of operator frequencies

    FAGeneticOperatorsMixture(
        operators::AbstractVector{<:GeneticOperator};
        c = 1E-2, eta = 1E-2, pmin = 0.01, pmax = 1.0
    ) = new(GeneticOperator[op for op in operators],
            FrequencyAdapter(length(operators);
                             c = c, eta = eta, pmin = pmin, pmax = pmax))
end

frequencies(opmix::FAGeneticOperatorsMixture) = frequencies(opmix.fa)

function next(opmix::FAGeneticOperatorsMixture)
    i = next(opmix.fa)
    return opmix.operators[i], i
end

function adjust!(opmix::FAGeneticOperatorsMixture, op_index::Int, candi_index::Int,
                 new_fitness::F, old_fitness::F, is_improved::Bool) where F
    # KLUDGE we don't know the fitness scheme, but if there is no improvement in fitness,
    # new_fitness==old_fitness, so it's ok to take the abs()
    update!(opmix.fa, op_index, is_improved ? abs(new_fitness-old_fitness) : 0.0)

    # also adjust the actual operator
    adjust!(opmix.operators[op_index], 0, candi_index, new_fitness, old_fitness, is_improved)
end

# FIXME use logging
function trace_state(io::IO, fa::FAGeneticOperatorsMixture, mode::Symbol)
    if mode == :verbose
        println(io, "FrequencyAdapting mixture rates:")
        for i in eachindex(fa.operators)
            @printf(io, " %s=%.2f", typeof(fa.operators[i]), fa.fa.p[i])
        end
    end
    println(io)
end
