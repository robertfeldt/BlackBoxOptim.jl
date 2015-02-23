abstract OptimizationProblem

type FixedDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}  # Objective functions
  ss::SearchSpace
  fmins::Union(Nothing, Vector{Float64})
  FixedDimProblem(name, funcs, ss, fmins = Nothing()) = begin
    new(name, funcs, ss, fmins)
  end
end

name(p::OptimizationProblem) = p.name

is_fixed_dimensional(p::Any) = false
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

fmin(p::OptimizationProblem) = fmins(p) != Nothing() ? fmins(p)[1] : nothing

ofunc(p::OptimizationProblem, index::Int) = p.funcs[index]

evalfunc(x, i::Int, p::OptimizationProblem) = ofunc(p, i)(x)

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

# Within ftol of a certain fmin
function fitness_is_within_ftol(p::OptimizationProblem, ftol, fitness; index = 1)
  fmins = BlackBoxOptim.fmins(p)
  (fmins == nothing) ? false : (abs(fmins[index] - fitness) < ftol)
end

type AnyDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}                 # Objective functions
  range_per_dimension::(Float64, Float64) # Default range per dimension
  fmins::Union(Nothing, Vector{Float64})
end

anydim_problem(name, f::Function, range, fmin::Float64) = AnyDimProblem(name, [f], range, [fmin])
anydim_problem(name, f::Function, range, fmins::Union(Nothing, Vector{Float64})) = AnyDimProblem(name, [f], range, fmins)
anydim_problem(name, f::Function, range) = AnyDimProblem(name, [f], range, nothing)

function as_fixed_dim_problem(p::AnyDimProblem, dim::Int)
  ss = symmetric_search_space(dim, p.range_per_dimension)
  FixedDimProblem(p.name, p.funcs, ss, p.fmins)
end

function as_fixed_dim_problem(p::FixedDimProblem, dim::Int)
  if numdims(p) != dim
    throw("Trying to set dimension $(dim) on a fixed dimensional problem of dimension $(numdims(p))")
  end
  p
end

function fixeddim_problem(f::Function; search_space = false, range = (-1.0, 1.0), 
  dims = 5, name = "unknown", fmins = Nothing())
  if search_space == false
    as_fixed_dim_problem(anydim_problem(name, f, range, fmins), dims)
  elseif typeof(search_space) <: SearchSpace
    FixedDimProblem(name, [f], search_space, fmins)
  elseif typeof(search_space) <: Array{(Float64,Float64),1}
    FixedDimProblem(name, [f], RangePerDimSearchSpace(search_space), fmins)
  else
    throw("Unknown search space $(search_space).")
  end
end

# A function set is specified through a dict mapping its function number
# to an optimization problem. We can create a fixed dimensional variant of
# an any dimensional function set with:
function as_fixed_dim_problem_set(ps::Dict{Any, Any}, dim::Int)
  as_fixed_dim_problem_set(ps, [dim])
end

# Create a fixed dim version of each problem in ps for each dim in dims.
function as_fixed_dim_problem_set(ps::Dict{Any, Any}, dims::Array{Int,1})
  next_free_index = 1
  result = Dict{Any, FixedDimProblem}()
  for(d in dims)
    for(p in values(ps))
      result[next_free_index] = as_fixed_dim_problem(p, d)
      next_free_index += 1
    end
  end
  result
end

include("single_objective.jl")