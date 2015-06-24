# Mutation clock operator is a more efficient way to mutate vectors than to generate
# a random value per variable in the vectors. It instead generates the number of variables
# to skip until the next mutation. Then it uses a sub-mutation operator to do the actual
# mutation. This is implemented as described in the paper:
#  Deb and Deb (2012), "Analyzing Mutation Schemes for Real-Parameter Genetic Algorithms"

num_vars_to_next_mutation_point(probMutation) = ceil( Int, (-log(rand())) / probMutation)

type MutationClock <: MutationOperator
  subMutationOperator::MutationOperator
  probMutation::Float64  # Probability of mutation of a variable
  nextVarToMutate::Int
  MutationClock(subMutationOperator, probMutation = 0.05) = new(subMutationOperator, 
    probMutation, num_vars_to_next_mutation_point(probMutation))
end

# Mutate each variable in a vector with the probability of mutation.
function apply{T <: Real}(mc::MutationClock, v::Vector{T})
  n = length(v)
  while mc.nextVarToMutate < n
    i = 1 + mc.nextVarToMutate
    v[i] = apply(mc.subMutationOperator, v[i])
    mc.nextVarToMutate += num_vars_to_next_mutation_point(mc.probMutation)
  end
  mc.nextVarToMutate -= n
  return v
end