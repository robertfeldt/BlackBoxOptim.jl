module Problems

export  OptimizationProblem,
        numdims, search_space

immutable OptimizationProblem
  name::ASCIIString

  # Objective functions
  funcs::Vector{Function}

  # True iff the problem can be instantiated in any number fo dimensions
  any_dimensional::Bool

  # Range per dimension
  range_per_dimension::(Float64, Float64)

  # Number of dimensions of the problem, if set.
  dimensions::Union(Nothing,Int)

  # Search space
  search_space::Union(Nothing,Array{(Float64, Float64)})
end

# Number of dimensions is by default 2, unless already specified.
numdims(p::OptimizationProblem) = (p.dimensions == nothing) ? 2 : p.dims

function search_space(p::OptimizationProblem)
  if (p.search_space == nothing)
    [p.range_per_dimension for i in 1:numdims(p)]
  else
    p.search_space
  end
end

function any_dimensional_problem(name, funcs, range_per_dimension)
  OptimizationProblem(name, funcs, true, range_per_dimension, nothing, nothing)
end

function one_dimensional_problem(name, funcs, range)
  OptimizationProblem(name, funcs, false, range, 1, [range])
end

examples = Dict{ASCIIString, OptimizationProblem}()

include("single_objective.jl")

end # module Problems