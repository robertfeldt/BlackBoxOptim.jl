# Mutation clock operator is a more efficient way to mutate vectors than to generate
# a random value per variable in the vectors. It instead generates the number of variables
# to skip until the next mutation. Then it uses a sub-mutation operator to do the actual
# mutation. This is based on the paper:
#  Deb and Deb (2012), "Analyzing Mutation Schemes for Real-Parameter Genetic Algorithms"
# but we use a Poisson distribution.

num_vars_to_next_mutation_point(probMutation) = ceil( Int, (-log(rand())) / probMutation)

# implements apply() that mutates one specific index of the parameter vector
abstract GibbsMutationOperator

# randomly mutate one index of parameter vector within its boundaries
immutable SimpleGibbsMutation <: GibbsMutationOperator
    ss::SearchSpace
end

search_space(m::SimpleGibbsMutation) = m.ss

# returns the mutated value of the specified dimension
function apply(mut::SimpleGibbsMutation, cur::Number, dim::Int)
  min, max = range_for_dim(mut.ss, dim)
  return min + rand() * (max-min)
end

type MutationClock{S<:GibbsMutationOperator} <: MutationOperator
  inner::S
  probMutation::Float64  # Probability of mutation of a variable
  nextVarToMutate::Int   # dimension index - 1
end

function MutationClock{S<:GibbsMutationOperator}(inner::S, probMutation::Float64 = 0.05)
  MutationClock{S}(inner, probMutation, num_vars_to_next_mutation_point(probMutation))
end

# Mutate each variable in a vector with the probability of mutation.
function apply!{T<:Real}(mc::MutationClock, v::Vector{T})
  n = length(v)
  while mc.nextVarToMutate < n
    i = 1 + mc.nextVarToMutate
    v[i] = apply(mc.inner, v[i], i)
    mc.nextVarToMutate += num_vars_to_next_mutation_point(mc.probMutation)
  end
  mc.nextVarToMutate -= n
  return v
end
