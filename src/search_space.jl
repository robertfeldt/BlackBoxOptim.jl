# A SearchSpace determines which candidate points can be
# considered in a search/optimization. This base class has very few restrictions
# and can allow varying number of dimensions etc.
abstract SearchSpace

# SearchSpace with a fixed number of dimensions. The vast majority of problems.
abstract FixedDimensionSearchSpace <: SearchSpace

# In a ContinuousSearchSpace each dimension has a continuous range of
# valid values.
abstract ContinuousSearchSpace <: FixedDimensionSearchSpace

# representative of the search space
typealias Individual Vector{Float64}

typealias ParamBounds @compat Tuple{Float64,Float64}

# Access individual range for a dim.
range_for_dim(css::ContinuousSearchSpace, i) = (mins(css)[i], maxs(css)[i])

ranges(css::ContinuousSearchSpace) = ParamBounds[range_for_dim(css, i) for i in 1:numdims(css)]

# Generate a number of individuals by random sampling in the search space.
function rand_individuals(css::ContinuousSearchSpace, numIndividuals)
  # Basically min + delta * rand(), individuals are stored in columns
  broadcast(+, mins(css), broadcast(*, deltas(css), rand(numdims(css), numIndividuals)))
end

# Generate a number of individuals via latin hypercube sampling. This should
# be the default when creating the initial population.
function rand_individuals_lhs(css::ContinuousSearchSpace, numIndividuals)
  Utils.latin_hypercube_sampling(mins(css), maxs(css), numIndividuals)
end

# Generate one random candidate.
function rand_individual(css::ContinuousSearchSpace)
  squeeze(rand_individuals(css, 1), 2)
end

# True iff ind is within the search space.
function Base.in(ind, css::ContinuousSearchSpace)
  @assert length(ind) == numdims(css)
  for i in eachindex(ind)
      if !(mins(css)[i] <= ind[i] <= maxs(css)[i])
        return false
      end
  end
  return true
end

# In a RangePerDimSearchSpace each dimension is specified as a range of valid
# values.
type RangePerDimSearchSpace <: ContinuousSearchSpace
  # We save the ranges as individual mins, maxs and deltas for faster access later.
  mins::Vector{Float64}
  maxs::Vector{Float64}
  deltas::Vector{Float64}

  function RangePerDimSearchSpace(ranges)
    mins = map(t -> t[1], ranges)
    maxs = map(t -> t[2], ranges)
    new(mins, maxs, (maxs - mins))
  end
end
mins(rss::RangePerDimSearchSpace) = rss.mins
maxs(rss::RangePerDimSearchSpace) = rss.maxs
deltas(rss::RangePerDimSearchSpace) = rss.deltas
numdims(rss::RangePerDimSearchSpace) = length(mins(rss))

diameters(rss::RangePerDimSearchSpace) = deltas(rss)

# Convenience function to create symmetric search spaces.
symmetric_search_space(numdims, range = (0.0, 1.0)) = RangePerDimSearchSpace(ParamBounds[range for i in 1:numdims])

# Create a feasible point (i.e. within the search space) given one which is
# outside.
feasible(v, ss::RangePerDimSearchSpace) = Float64[ clamp( v[i], mins(ss)[i], maxs(ss)[i] ) for i in eachindex(v) ]
