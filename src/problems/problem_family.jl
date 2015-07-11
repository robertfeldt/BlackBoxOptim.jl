# family of FunctionBasedProblem parameterized by the search space dimension
type FunctionBasedProblemFamily{F,FS<:FitnessScheme}
  objfunc::Function                          # Objective function
  name::ASCIIString
  fitness_scheme::FS
  range_per_dimension::ParamBounds        # Default range per dimension
  opt_value::Nullable{F}                  # optional optimal value

  FunctionBasedProblemFamily(objfunc::Function, name::ASCIIString, fitness_scheme::FS, range::ParamBounds,
                opt_value::Nullable{F}) = new(objfunc, name, fitness_scheme, range, opt_value)
end

Base.convert{F,FS<:FitnessScheme}(FunctionBasedProblemFamily, objfunc::Function, name::ASCIIString, fitness_scheme::FS, range::ParamBounds,
                opt_value::Nullable{F}) = FunctionBasedProblemFamily{F,FS}(objfunc, name, fitness_scheme, range, opt_value)

objfunc(family::FunctionBasedProblemFamily) = family.objfunc

# generates the problem with the specified number of dimensions
function fixed_dim_problem(family::FunctionBasedProblemFamily, ndim::Int)
  ss = symmetric_search_space(ndim, family.range_per_dimension)
  if isnull(family.opt_value)
    return convert(FunctionBasedProblem, family.objfunc, family.name, family.fitness_scheme, ss)
  else
    return convert(FunctionBasedProblem, family.objfunc, family.name, family.fitness_scheme, ss, get(family.opt_value))
  end
end

function MinimizationProblemFamily(f::Function, name::ASCIIString, range::ParamBounds, fmin::Float64)
  convert(FunctionBasedProblemFamily, f, name, ScalarFitness{true}(), range, convert(Nullable{Float64}, fmin))
end

function MinimizationProblemFamily(f::Function, name::ASCIIString, range::ParamBounds)
  convert(FunctionBasedProblemFamily, f, name, ScalarFitness{true}(), range, Nullable{Float64}())
end

function minimization_problem(f::Function, name::ASCIIString, range::ParamBounds, ndim::Int)
    fixed_dim_problem(MinimizationProblemFamily(f, name, range), ndim)
end

function minimization_problem(f::Function, range::ParamBounds, ndim::Int)
    fixed_dim_problem(MinimizationProblemFamily(f, "unknown name", range), ndim)
end

function minimization_problem(f::Function, name::ASCIIString, range::ParamBounds, ndim::Int, fmin::Float64)
    fixed_dim_problem(MinimizationProblemFamily(f, name, range, fmin), ndim)
end

# A function set is specified through a dict mapping its function number
# to an optimization problem. We can create a fixed dimensional variant of
# an any dimensional function set with:
function problem_set(ps::Dict{Any, Any}, dim::Int)
  problem_set(ps, [dim])
end

# Create a fixed dim version of each problem in ps for each dim in dims.
function problem_set(ps::Dict{Any, FunctionBasedProblemFamily}, dims::Vector{Int})
  next_free_index = 1
  result = Dict{Int, FunctionBasedProblem}()
  for d in dims
    for p in values(ps)
      result[next_free_index] = generate(p, d)
      next_free_index += 1
    end
  end
  result
end
