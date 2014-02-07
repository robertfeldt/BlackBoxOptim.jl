# This is the R2-IBEA weight vector generation algorithm as described in:
# Dung H. Phan and Junichi Suzuki, "R2-IBEA: R2 Indicator Based Evolutionary
# Algorithm for Multiobjective Optimization", 2013.

function generate_random_weight_vector(m)

  x = sort(rand(m-1))

  w = zeros(m)
  w[1] = x[1]
  for i in 2:(m-1)
    w[i] = x[i] - x[i-1]
  end
  w[m] = 1 - x[m-1]

  w

end

function generate_weight_vectors(numobjectives, numvectors, maxnum_iterations = 10000)

  # The final returned set of vectors will be returned as a 
  # numobjectives * numvectors matrix but we need room for one more column
  # to keep the last candidate vector to be added.
  generated_vectors = zeros(numobjectives, numvectors+1) # W
  numgenerated = 0

  for t in 1:maxnum_iterations

    numgenerated += 1
    generated_vectors[:, numgenerated] = generate_random_weight_vector(numobjectives)

    if numgenerated > numvectors

      i = index_to_column_with_min_hypervolume_contribution(generated_vectors)
      generated_vectors[:, i] = generated_vectors[:, numgenerated]
      numgenerated -= 1

    end
  end

  # Delete the additional column and return the correct number of vectors
  generated_vectors[:, 1:numobjectives]

end

function index_to_column_with_min_hypervolume_contribution(vectors)

end


#
# Calculate the approximate hypervolume contribution of each of the L
# M-dimensional fitness vectors given as a M*L matrix, v. The calculation is
# via monte-carlo simulation.
#
function approximate_hypervolume_monte_carlo(v::Array{Float64, 2}; 
  upperbound = maximum(v, 2), numsamples = 100000)

  M, L = size(v)

  # Get the lower bound point and the delta, i.e. size of box between lower and upper bounds.
  lowerbound = minimum(v, 2)
  delta = upperbound - lowerbound

  # Alloc mem once so we can write into it later for speed.
  sampledpoint = zeros(M)

  # Count the number of sampled points that are dominated by the given points.
  num_dominated = 0
  for i in 1:numsamples

    # Sample in the box from lowerbound to upperbound.
    for m in 1:M
      sampledpoint[m] = lowerbound[m] + rand() * delta[m]
    end

    # Compare sampled point to points on the front and count if it is dominated.
    for frontpoint in 1:L
      if dominates(v[:, frontpoint], sampledpoint)
        num_dominated += 1
        break
      end
    end

  end

  # Approximate hypervolume based on how many of the sampled points were dominated.
  return prod(delta) * num_dominated / numsamples

end

#
# Return true iff points a dominates point b in the Pareto sense.
#
function dominates(a, b)

  dominates_in_at_least_one = false

  for i in 1:length(a)
    if a[i] > b[i]
      return false
    elseif a[i] < b[i]
      dominates_in_at_least_one = true
    end
  end

  return dominates_in_at_least_one

end
