abstract OptimizationProblem

type FixedDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}  # Objective functions
  ss::SearchSpace
  fmins::Union(Nothing, Vector{Float64})
end

is_fixed_dimensional(p::OptimizationProblem) = false
is_fixed_dimensional(p::FixedDimProblem) = true

is_any_dimensional(p::OptimizationProblem) = !is_fixed_dimensional(p)

numfuncs(p::OptimizationProblem) = length(p.funcs)
is_single_objective_problem(p::OptimizationProblem) = numfuncs(p) == 1

is_multi_objective_problem(p::OptimizationProblem) = !is_single_objective_problem(p)

numdims(p::OptimizationProblem) = nothing
numdims(p::FixedDimProblem) = numdims(p.ss)

search_space(p::OptimizationProblem) = nothing
search_space(p::FixedDimProblem) = p.ss

fmins(p::OptimizationProblem) = p.fmins

ofunc(p::OptimizationProblem, index::Int64) = p.funcs[index]

evalfunc(x, i::Int64, p::OptimizationProblem) = ofunc(p, i)(x)

# Evaluate fitness of a candidate solution on the 1st objective function of a problem.
eval1(x, p::OptimizationProblem) = evalfunc(x, 1, p)

# Evaluate fitness of a candidate solution on all objective functions of a problem.
evalall(x, p::OptimizationProblem) = begin
  n = length(p.funcs)
  results = zeros(n)
  for(i in 1:n)
    results[i] = evalfunc(x, i, p)
  end
  results
end

type AnyDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}                 # Objective functions
  range_per_dimension::(Float64, Float64) # Default range per dimension
  fmins::Union(Nothing, Vector{Float64})
end

anydim_problem(name, f::Function, range, fmin::Float64) = AnyDimProblem(name, [f], range, [fmin])
anydim_problem(name, f::Function, range) = AnyDimProblem(name, [f], range, nothing)

function as_fixed_dim_problem(p::AnyDimProblem, dim::Int64)
  ss = symmetric_search_space(dim, p.range_per_dimension)
  FixedDimProblem(p.name, p.funcs, ss, p.fmins)
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