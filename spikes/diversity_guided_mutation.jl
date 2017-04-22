# This implements the DGM operator for Diversity-Guided Mutation as described in
#  
#  M. S. Alam, Md. M. Islam, X. Yao, K. Murase, "Diversity Guided Evolutionary 
#  Programming: A novel approach for continuous optimization", ASOC, 2013
#

# DGM operator - this is a slight variation from the one in their paper
# were we do not need to set a min distance but rather select the min and non-min
# sets by randomly sampling a subset of individuals in the population.
#
# Input:
#  - population is a matrix of Float64 where each column represents one individual in the population
#  - parent is an Int index that points out one column in the population which will act as the parent to be mutated
#  - population_size
# Output:
#  - a child that is a diversity-mutated variant of the parent
#
function diversity_guided_mutation(population::Array{Float64, 2}, parent::Int; 
  sample_size = 10,
  indices = collect(1:size(population, 2)))

  p = population[:,parent]
  population_size = size(population, 2)
  shuffle!(indices)

  # Calculate the distances from sample_size randomly sampled individuals
  # to the parent.
  distances = zeros(sample_size)
  subset = zeros(Int, sample_size)
  di = i = 1
  while di <= sample_size
    si = indices[i]
    if si == parent
      i += 1
      si = indices[i]
    end
    distances[di] = norm(population[:,si] - p)
    di += 1
    i += 1
  end
  
  # Sort the distances in the subset so we can have the low 3rd as Nx and the 
  # high 3rd as Pop-Nx.
  perm = sortperm(distances)
  third_size = round(sample_size/3)
  subset = subset[perm]
  nx = subset[1:third_size]
  pt_minus_nx = subset[]