"""
  The base abstract type for all optimization problems.
  `FS` is a type of a problem's `FitnessScheme`.
"""
abstract OptimizationProblem{FS<:FitnessScheme}

# common definitions for `OptimizationProblem`
# (enforce field names of subtypes)
name(p::OptimizationProblem) = p.name
fitness_scheme(p::OptimizationProblem) = p.fitness_scheme
fitness_type(p::OptimizationProblem) = fitness_type(fitness_scheme(p))
numobjectives(p::OptimizationProblem) = numobjectives(fitness_scheme(p))
search_space(p::OptimizationProblem) = p.ss
numdims(p::OptimizationProblem) = numdims(search_space(p))

"""
  Checks if the current best fitness is within the given range
  of a known best fitness.
  By default the best fitness is unknown.
"""
fitness_is_within_ftol(p::OptimizationProblem, fitness, atol::Number) = false

"""
  `OptimizationProblem` with the objective function
  defined by some Julia `Function`.
"""
abstract FunctionBasedProblem{FS<:FitnessScheme} <: OptimizationProblem{FS}

objfunc(p::FunctionBasedProblem) = p.objfunc

"""
  Evaluates fitness of a candidate.
"""
fitness(x, p::FunctionBasedProblem) = objfunc(p)(x)

"""
  Problem with unknown global optimum with fitness
  defined by some Julia `Function`.
"""
type UnboundedProblem{FS<:FitnessScheme, SS<:SearchSpace} <: FunctionBasedProblem{FS}
  objfunc::Function  # Objective function
  name::ASCIIString
  fitness_scheme::FS
  ss::SS          # search space

  function UnboundedProblem(objfunc::Function, name::ASCIIString, fitness_scheme::FS, ss::SS)
    new(objfunc, name, fitness_scheme, ss)
  end
end

Base.copy{FS,SS}(p::UnboundedProblem{FS,SS}) =
  UnboundedProblem{FS,SS}(p.objfunc, copy(p.name), p.fitness_scheme, p.ss)

"""
  Problem with known global optimum,
  defined by some Julia `Function`.
"""
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

Base.copy{F,FS,SS}(p::BoundedProblem{F,FS,SS}) =
  BoundedProblem{F,FS,SS}(p.objfunc, copy(p.name), p.fitness_scheme, p.ss, p.opt_value)

opt_value(p::BoundedProblem) = p.opt_value
fitness_is_within_ftol(p::BoundedProblem, fitness, atol::Number) = norm(opt_value(p) - fitness) < atol

function Base.convert{FS<:FitnessScheme,SS<:SearchSpace}(::Type{FunctionBasedProblem}, func::Function, name::ASCIIString, fitness_scheme::FS, ss::SS)
  return UnboundedProblem{FS,SS}(func, name, fitness_scheme, ss)
end

function Base.convert{F,FS<:FitnessScheme,SS<:SearchSpace}(::Type{FunctionBasedProblem}, func::Function, name::ASCIIString, fitness_scheme::FS, ss::SS, opt_value::F)
  return BoundedProblem{F,FS,SS}(func, name, fitness_scheme, ss, opt_value)
end
