num_vars_to_next_mutation_point(probMutation) = ceil( Int, (-log(rand())) / probMutation)

"""
    Provides `apply()` operator that mutates one specified dimension of a parameter
    vector.
"""
abstract GibbsMutationOperator <: MutationOperator

# apply operator to each dimension
function apply!{T<:Real}(m::GibbsMutationOperator, params::AbstractVector{T}, target_index::Int)
    @inbounds @simd for i in eachindex(params)
        params[i] = apply(m, params[i], i, target_index)
    end
    return params
end

"""
    Uniform mutation of a parameter vector.
"""
immutable UniformMutation{SS<:SearchSpace} <: GibbsMutationOperator
    ss::SS

    @compat (::Type{UniformMutation}){SS<:SearchSpace}(ss::SS) = new{SS}(ss)
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
type MutationClock{S<:GibbsMutationOperator} <: MutationOperator
    inner::S
    rate::Float64           # Probability to mutate a dimension
    nextVarToMutate::Int    # dimension index - 1

    @compat (::Type{MutationClock}){S<:GibbsMutationOperator}(inner::S, rate::Float64 = 0.05) =
        new{S}(inner, rate, 1 + num_vars_to_next_mutation_point(rate))
end

function apply!{T<:Real}(mc::MutationClock, v::AbstractVector{T}, target_index::Int)
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
