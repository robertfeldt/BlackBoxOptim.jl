"""
Several genetic operators chained together.
"""
struct OperatorPipeline <: GeneticOperator
    ops::Vector{GeneticOperator}

    function OperatorPipeline(operators...)
        ensure_nargs_match(operators) || throw("The number of ")
        new(operators)
    end
end

function apply(p::OperatorPipeline, parents::AbstractVector{<:AbstractVector{<:Real}})
    for i in 1:length(p.ops)
        parents = apply(p.ops[i], parents)
    end
    return parents
end

match(o1::GeneticOperator, o2::GeneticOperator) = numchildren(o1) == numparents(o2)
# Mutation ops can accept any number of parents
match(o1::GeneticOperator, o2::MutationOperator) = true

function ensure_nargs_match(operators)
    for i in 1:(length(operators)-1)
        match(operators[i], operators[i+1]) || return false
    end
    return true
end
