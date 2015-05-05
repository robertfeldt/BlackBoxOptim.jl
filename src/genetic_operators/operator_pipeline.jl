# Several genetic operators can be chained together.

type OperatorPipeline <: GeneticOperator
  ops::Vector{GeneticOperator}
  OperatorPipeline(operators...) = begin
    if !ensure_nargs_match(operators)
      throw("The number of ")
    end
    new(operators)
  end
end

function apply{T <: Real}(p::OperatorPipeline, parents::Vector{Vector{T}})
  for i in 1:length(p.ops)
    parents = apply(p.ops[i], parents)
  end
  parents
end

match(o1::GeneticOperator, o2::GeneticOperator) = numchildren(o1) == numparents(o2)
# Mutation ops can accept any number of parents
match(o1::GeneticOperator, o2::MutationOperator) = true

function ensure_nargs_match(operators)
  for i in 1:(length(operators)-1)
    if !match(operators[i], operators[i+1])
      return false
    end
  end
  true
end