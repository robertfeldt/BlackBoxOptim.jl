using Distributions

CuckooSearchDefaultOptions = {
  "AlienEggsDiscoveryRate" => 0.25,
}

# This is based on the "Engineering Optimisation by Cuckoo Search" paper
# by Yang and Deb 2012. IMHO, the paper lacks details and it is hard to know
# which parts of what is described are essential and which is accidental.
# Below we try to stay quite close to the original paper. However, what is 
# called nests in the cuckoo search papers by Yang seems to be basically a 
# population of candidate solutions. Thus, we just call it a population.
type CuckooSearch <: PopulationOptimizer
  name::ASCIIString

  # A population is a matrix of floats.
  population::Array{Float64, 2}

  # A search space is defined by the min and max values (in tuples) for each
  # of its dimenions. The dimension is the length of an individual, i.e. the
  # number of Float64 values in it.
  search_space::SearchSpace

  # Options
  options::Dict{Any,Any}

  function CuckooSearch(name, pop, ss, options)
    new(name, pop, ss, merge(CuckooSearchDefaultOptions, options))
  end
end

