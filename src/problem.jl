"""
The base abstract type for all optimization problems.
`FS` is a type of a problem's `FitnessScheme`.
"""
abstract type OptimizationProblem{FS<:FitnessScheme} end

# common definitions for `OptimizationProblem`
# (enforce field names of subtypes)
name(p::OptimizationProblem) = p.name
fitness_scheme_type(::Type{<:OptimizationProblem{FS}}) where FS = FS
fitness_type(::Type{P}) where P <: OptimizationProblem = fitness_type(fitness_scheme_type(P))
fitness_scheme(p::OptimizationProblem) = p.fitness_scheme
fitness_type(p::P) where P <: OptimizationProblem = fitness_type(P)
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
defined by some Julia `Function` and search space.

Optionally, a known optimal value could be provided to terminate
the optimization once it is reached.
"""
mutable struct FunctionBasedProblem{FS<:FitnessScheme,SS<:SearchSpace,FO} <: OptimizationProblem{FS}
    objfunc::Function     # Objective function
    name::String
    fitness_scheme::FS
    ss::SS                # search space
    opt_value::FO         # known optimal value or nothing

    function FunctionBasedProblem(objfunc::Function, name::String,
                                  fitness_scheme::FS, ss::SS,
                                  opt_value::FO = nothing) where {FS<:FitnessScheme,SS<:SearchSpace,FO}
        if FO <: Number
            fitness_type(fitness_scheme) == FO ||
                throw(ArgumentError("Fitness type ($(fitness_type(fitness_scheme))) and opt_value type ($(FO)) do not match"))
            if isnafitness(opt_value, fitness_scheme) # if NA fitness is given, the problem is unbouned
                opt_value = nothing
                #throw(ArgumentError("Known optimal value cannot be NA"))
            end
        end
        new{FS,SS,typeof(opt_value)}(objfunc, name, fitness_scheme, ss, opt_value)
    end
end

objfunc(p::FunctionBasedProblem) = p.objfunc

"""
    fitness(x, p::OptimizationProblem)

Evaluates the fitness of a candidate.
"""
fitness(x, p::FunctionBasedProblem) = objfunc(p)(x)

Base.copy(p::FunctionBasedProblem) =
    FunctionBasedProblem(deepcopy(p.objfunc), deepcopy(p.name),
                         p.fitness_scheme, p.ss, p.opt_value)

opt_value(p::FunctionBasedProblem) = p.opt_value

# checks if the optimal fitness is reached for bounded scalar problem
fitness_is_within_ftol(p::FunctionBasedProblem{<:FitnessScheme, <:SearchSpace, <:Number},
                       fitness, atol::Number) = norm(opt_value(p) - fitness) < atol

# no check for unbounded problem (by default, see the definition for OptimizationProblem above)
# fitness_is_within_ftol{FS,SS}(p::FunctionBasedProblem{FS,SS,Void}, fitness::Number, atol::Number) = true
