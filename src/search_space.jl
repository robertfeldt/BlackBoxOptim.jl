# A SearchSpace determines which candidate points can be 
# considered in a search/optimization. This base class has very few restrictions
# and can allow varying number of dimensions etc.
abstract SearchSpace

# SearchSpace with a fixed number of dimensions. The vast majority of problems.
abstract FixedDimensionSearchSpace <: SearchSpace

# In a ContinuousSearchSpace each dimension has a continuous range of
# valid values.
abstract ContinuousSearchSpace <: FixedDimensionSearchSpace

numdims(css::ContinuousSearchSpace) = css.numdims
ranges(css::ContinuousSearchSpace) = [range_for_dim(css, i) for i in 1:numdims(css)]

# Access individual range for a dim.
range_for_dim(css::ContinuousSearchSpace, i::Int) = (mins(css)[i], maxs(css)[i])

# Generate a number of individuals by random sampling in the search space.
function rand_individuals(css::ContinuousSearchSpace, numIndividuals)
  # Basically min + delta * rand(), but broadcast over the columns...
  broadcast(+, mins(css), broadcast(*, deltas(css), rand(numIndividuals, numdims(css))))
end

# Generate a number of individuals via latin hypercube sampling. This should
# be the default when creating the initial population.
function rand_individuals_lhs(css::ContinuousSearchSpace, numIndividuals)
  BlackBoxOptim.Utils.latin_hypercube_sampling(mins(css), maxs(css), numIndividuals)
end

# Generate one random candidate.
function rand_individual(css::ContinuousSearchSpace)
  rand_individuals(css, 1)[1,:]
end

# True iff ind is within the search space.
function isinspace(ind, css::ContinuousSearchSpace)
  all(mins(css) .<= ind .<= maxs(css))
end

# In a RangePerDimSearchSpace each dimension is specified as a range of valid
# values.
type RangePerDimSearchSpace <: ContinuousSearchSpace
  numdims::Int
  # We save the ranges as individual mins, maxs and deltas for faster access later.
  mins::Array{Float64,2}
  maxs::Array{Float64,2}
  deltas::Array{Float64,2}

  function RangePerDimSearchSpace(ranges)
    mins = map(t -> t[1], ranges)'
    maxs = map(t -> t[2], ranges)'
    new(length(ranges), mins, maxs, (maxs - mins))
  end
end
mins(rss::RangePerDimSearchSpace) = rss.mins
maxs(rss::RangePerDimSearchSpace) = rss.maxs
deltas(rss::RangePerDimSearchSpace) = rss.deltas

#convert(::Type{ContinuousSearchSpace}, ranges::Array{(Float64,Float64),1}) =
#  RangePerDimSearchSpace(ranges)

# Convenience function to create symmetric search spaces.
symmetric_search_space(numdims, range = (0.0, 1.0)) = RangePerDimSearchSpace([range for i in 1:numdims])

# Create a feasible point (i.e. within the search space) given one which is
# outside.
function feasible(v, ss::RangePerDimSearchSpace)
  minimum(hcat(maxs(ss)', maximum(hcat(mins(ss)', v), 2)), 2)
end
