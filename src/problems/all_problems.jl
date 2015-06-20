abstract OptimizationProblem

type FixedDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}  # Objective functions
  ss::SearchSpace
  fmins::Nullable{Vector{Float64}}

  function FixedDimProblem(name, funcs, ss, fmins::Nullable{Vector{Float64}} = Nullable{Vector{Float64}}())
    @assert isnull(fmins) || length(funcs) == length(get(fmins))
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

fmin(p::OptimizationProblem) = !isnull(fmins(p)) ? Nullable{Float64}( get(fmins(p))[1] ) : Nullable{Float64}()

ofunc(p::OptimizationProblem, i) = p.funcs[i]

evalfunc(x, i, p::OptimizationProblem) = ofunc(p, i)(x)

# Evaluate fitness of a candidate solution on the 1st objective function of a problem.
eval1(x, p::OptimizationProblem) = evalfunc(x, 1, p)

# Evaluate fitness of a candidate solution on all objective functions of a problem.
evalall(x, p::OptimizationProblem) = Float64[ f(x) for f in p.funcs ]

# Within ftol of a certain fmin
fitness_is_within_ftol(p::OptimizationProblem, ftol, fitness, index = 1) = isnull(fmins(p)) ? false : ( abs(get(fmins(p))[index] - fitness) < ftol )

type AnyDimProblem <: OptimizationProblem
  name::ASCIIString
  funcs::Vector{Function}                 # Objective functions
  range_per_dimension::ParamBounds        # Default range per dimension
  fmins::Nullable{Vector{Float64}}
end

function Base.convert( ::Type{Nullable{Vector{Float64}}}, v::Nullable{Float64} )
  isnull(v) ? Nullable{Vector{Float64}}() : Nullable{Vector{Float64}}( Float64[get(v)] )
end

function Base.convert( ::Type{Nullable{Vector{Float64}}}, v::Float64 )
  Nullable{Vector{Float64}}( Float64[v] )
end

anydim_problem(name, f::Function, range, fmin = Float64 ) = AnyDimProblem( name, Function[f], range, convert( Nullable{Vector{Float64}}, fmin ) )
anydim_problem(name, f::Function, range ) = AnyDimProblem( name, Function[f], range, Nullable{Vector{Float64}}() )

function as_fixed_dim_problem(p::AnyDimProblem, dim::Int)
  ss = symmetric_search_space(dim, p.range_per_dimension)
  FixedDimProblem(p.name, p.funcs, ss, p.fmins)
end

function as_fixed_dim_problem(p::FixedDimProblem, dim::Int)
  if numdims(p) != dim
    throw(DimensionMismatch("Trying to set dimension $(dim) on a fixed dimensional problem of dimension $(numdims(p))"))
  end
  p
end

function fixeddim_problem(f::Function; search_space = false, range = (-1.0, 1.0),
  dims = 5, name = "unknown", fmin = Nullable{Float64}() )
  fmins = convert( Nullable{Vector{Float64}}, fmin )
  if search_space == false
    as_fixed_dim_problem(anydim_problem(name, f, range, fmins), dims)
  elseif typeof(search_space) <: SearchSpace
    FixedDimProblem(name, Function[f], search_space, fmins)
  elseif typeof(search_space) <: Vector{ParamBounds}
    FixedDimProblem(name, Function[f], RangePerDimSearchSpace(search_space), fmins)
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
