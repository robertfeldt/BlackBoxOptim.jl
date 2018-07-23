"""
Wraps the mutation operator, so that it could
be used as crossover operator.
"""
struct MutationWrapper{OP<:MutationOperator} <: CrossoverOperator{1,1}
    inner::OP

    MutationWrapper(mutation::OP) where {OP<:MutationOperator} =
        new{OP}(mutation)
end

function apply!(wrapper::MutationWrapper,
                target::Individual, targetIndex::Int, pop, parentIndices)
    @assert length(parentIndices) == 1
    apply!(wrapper.inner, copyto!(target, viewer(pop, parentIndices[1])), targetIndex)
end

Base.show(io::IO, xover::MutationWrapper) = print(io, "MutationWrapper(", xover.inner, ")")

# convert mutation to crossover by wrapping it in `MutationWrapper`
Base.convert(::Type{CrossoverOperator}, op::MutationOperator) = MutationWrapper(op)
