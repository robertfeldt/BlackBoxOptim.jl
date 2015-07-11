# root type for all optimization problems
abstract OptimizationProblem{FS<:FitnessScheme}

# common definitions for OptimizationProblem
# (enforce field names of subtypes)
name(p::OptimizationProblem) = p.name
fitness_scheme(p::OptimizationProblem) = p.fitness_scheme
fitness_type(p::OptimizationProblem) = fitness_type(fitness_scheme(p))
numobjectives(p::OptimizationProblem) = numobjectives(fitness_scheme(p))
search_space(p::OptimizationProblem) = p.ss
numdims(p::OptimizationProblem) = numdims(search_space(p))

# by default the optimum is unknown
fitness_is_within_ftol(p::OptimizationProblem, fitness, atol::Number) = false

# problem with the objective function defined by some Julia function
abstract FunctionBasedProblem{FS<:FitnessScheme} <: OptimizationProblem{FS}

objfunc(p::FunctionBasedProblem) = p.objfunc

# Evaluate fitness of a candidate
fitness(x, p::FunctionBasedProblem) = objfunc(p)(x)

# problem with unknown global optimum,
# a wrapper around Julia objective function
type UnboundedProblem{FS<:FitnessScheme, SS<:SearchSpace} <: FunctionBasedProblem{FS}
  objfunc::Function  # Objective function
  name::ASCIIString
  fitness_scheme::FS
  ss::SS          # search space

  function UnboundedProblem(objfunc::Function, name::ASCIIString, fitness_scheme::FS, ss::SS)
    new(objfunc, name, fitness_scheme, ss)
  end
end

# problem with known global optimum,
# a wrapper around Julia objective function
type BoundedProblem{F, FS<:FitnessScheme, SS<:SearchSpace} <: FunctionBasedProblem{FS}
  objfunc::Function  # Objective function
  name::ASCIIString
  fitness_scheme::FS
  ss::SS          # search space
  opt_value::F    # known optimal value

  function BoundedProblem(objfunc::Function, name::ASCIIString, fitness_scheme::FS, ss::SS, opt_value::F)
    if isnafitness(opt_value, fitness_scheme)
      throw(ArgumentError("The object function must be specified"))
    end
    new(objfunc, name, fitness_scheme, ss, opt_value)
  end
end

opt_value(p::BoundedProblem) = p.opt_value
fitness_is_within_ftol(p::BoundedProblem, fitness, atol::Number) = norm(opt_value(p) - fitness) < atol

function Base.convert{FS<:FitnessScheme,SS<:SearchSpace}(::Type{FunctionBasedProblem}, func::Function, name::ASCIIString, fitness_scheme::FS, ss::SS)
  return UnboundedProblem{FS,SS}(func, name, fitness_scheme, ss)
end

function Base.convert{F,FS<:FitnessScheme,SS<:SearchSpace}(::Type{FunctionBasedProblem}, func::Function, name::ASCIIString, fitness_scheme::FS, ss::SS, opt_value::F)
  return BoundedProblem{F,FS,SS}(func, name, fitness_scheme, ss, opt_value)
end
