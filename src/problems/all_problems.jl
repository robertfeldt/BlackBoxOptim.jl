module Problems

import BlackBoxOptim.numdims, # Since we will add to it
        BlackBoxOptim.SearchSpace, BlackBoxOptim.symmetric_search_space

export  OptimizationProblem,
        numdims, search_space, set_numdims!

abstract OptimizationProblem

type FixedDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}                 # Objective functions
  ss::SearchSpace  
end

is_fixed_dimensional(p::OptimizationProblem) = false
is_fixed_dimensional(p::FixedDimProblem) = true

is_any_dimensional(p::OptimizationProblem) = not(is_fixed_dimensional(p))

is_single_objective_problem(p::OptimizationProblem) = length(p.funcs) == 1

is_multi_objective_problem(p::OptimizationProblem) = not(is_single_objective_problem(p))

numdims(p::OptimizationProblem) = nothing
numdims(p::FixedDimProblem) = numdims(p.ss)

search_space(p::OptimizationProblem) = nothing
search_space(p::FixedDimProblem) = p.ss

type AnyDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}                 # Objective functions
  range_per_dimension::(Float64, Float64) # Default range per dimension 
end

anydim_problem(name, f::Function, range) = AnyDimensionalProblem(name, [f], range)

function as_fixed_dim_problem(p::AnyDimProblem, dim::Int64)
  ss = symmetric_search_space(dim, p.range_per_dimension)
  FixedDimProblem(p.name, p.funcs, ss)
end

# A function set is specified through a duct mapping its function number
# to an optimization problem. We can create a fixed dimensional variant of
# an any dimensional function set with:
function as_fixed_dim_problem_set(ps::Dict{Any, AnyDimProblem}, dim::Int64)
  result = Dict{Any, FixedDimProblem}()
  for(key in keys(ps))
    result[key] = as_fixed_dim_problem(ps[key], dim)
  end
  result
end

include("single_objective.jl")

end # module Problems