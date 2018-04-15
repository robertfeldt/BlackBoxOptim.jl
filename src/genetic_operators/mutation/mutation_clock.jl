num_vars_to_next_mutation_point(probMutation) = ceil( Int, (-log(rand())) / probMutation)

"""
Provides `apply()` operator that mutates one specified dimension of a parameter
vector.
"""
abstract type GibbsMutationOperator <: MutationOperator end

# apply operator to each dimension
function apply!(m::GibbsMutationOperator, params::AbstractVector{<:Real}, target_index::Int)
    @inbounds @simd for i in eachindex(params)
        params[i] = apply(m, params[i], i, target_index)
    end
    return params
end

"""
Uniform mutation of a parameter vector.
"""
struct UniformMutation{SS<:SearchSpace} <: GibbsMutationOperator
    ss::SS

    UniformMutation(ss::SS) where {SS<:SearchSpace} = new{SS}(ss)
end

search_space(m::UniformMutation) = m.ss

@inline apply(m::UniformMutation, v::Number, dim::Int, target_index::Int) =
    return (mins(m.ss)[dim] + rand() * deltas(m.ss)[dim])

"""
Mutation clock operator is a more efficient way to mutate vectors than to generate
a random value per variable in the vectors. It instead generates the number of variables
to skip until the next mutation. Then it uses a sub-mutation operator to do the actual
mutation. This is based on the paper:
    Deb and Deb (2012), "Analyzing Mutation Schemes for Real-Parameter Genetic Algorithms"
but we use a Poisson distribution.
"""
mutable struct MutationClock{S<:GibbsMutationOperator} <: MutationOperator
    inner::S
    rate::Float64           # Probability to mutate a dimension
    nextVarToMutate::Int    # dimension index - 1

    MutationClock(inner::S, rate::Float64 = 0.05) where {S<:GibbsMutationOperator} =
        new{S}(inner, rate, 1 + num_vars_to_next_mutation_point(rate))
end

function apply!(mc::MutationClock, v::AbstractVector{<:Real}, target_index::Int)
    n = length(v)
    i = mc.nextVarToMutate
    @inbounds while i <= n
        v[i] = apply(mc.inner, v[i], i, target_index)
        i += num_vars_to_next_mutation_point(mc.rate)
    end
    @assert i > n
    mc.nextVarToMutate = i - n
    return v
end
