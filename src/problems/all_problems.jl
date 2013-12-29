abstract OptimizationProblem

type FixedDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}  # Objective functions
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

# Evaluate fitness of a candidate solution on the 1st objective function of a problem.
eval1(x, p::OptimizationProblem) = p.funcs[1](x)

type AnyDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}                 # Objective functions
  range_per_dimension::(Float64, Float64) # Default range per dimension 
end

anydim_problem(name, f::Function, range) = AnyDimProblem(name, [f], range)

function as_fixed_dim_problem(p::AnyDimProblem, dim::Int64)
  ss = symmetric_search_space(dim, p.range_per_dimension)
  FixedDimProblem(p.name, p.funcs, ss)
end

# A function set is specified through a duct mapping its function number
# to an optimization problem. We can create a fixed dimensional variant of
# an any dimensional function set with:
function as_fixed_dim_problem_set(ps::Dict{Any, Any}, dim::Int64)
  result = Dict{Any, FixedDimProblem}()
  for(key in keys(ps))
    result[key] = as_fixed_dim_problem(ps[key], dim)
  end
  result
end

include("single_objective.jl")