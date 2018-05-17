
"""
Base class for problem families.

It is an abstraction for problem parameterization (e.g by the
number of the search space dimensions) that allows
to instantiate `OptimizationProblem` for the concrete parameters.
"""
abstract type ProblemFamily{FS<:FitnessScheme} end

"""
Family of `FunctionBasedProblem` optimization problems
parameterized by the number of search space dimensions.
"""
struct FunctionBasedProblemFamily{F,FS<:FitnessScheme,FO} <: ProblemFamily{FS}
    objfunc::Function                     # Objective function
    name::String
    fitness_scheme::FS
    reserved_ss::RangePerDimSearchSpace   # search space for the first reserved dimensions
    range_per_dim::ParamBounds            # Default range per dimension
    opt_value::FO                         # optional optimal value, or nothing

    function FunctionBasedProblemFamily(
            objfunc::Function, name::String,
            fitness_scheme::FS, range::ParamBounds, opt_value::FO = nothing,
            reserved_ss::RangePerDimSearchSpace = ZERO_SEARCH_SPACE
    ) where {FS<:FitnessScheme, FO}
        if FO <: Number
            fitness_type(fitness_scheme) == FO ||
                throw(ArgumentError("Fitness type ($(fitness_type(fitness_scheme))) and opt_value type ($(FO)) do not match"))
        end
        new{FO, FS, FO}(objfunc, name, fitness_scheme,
                        reserved_ss, range, opt_value)
    end
end

objfunc(family::FunctionBasedProblemFamily) = family.objfunc

"""
    instantiate(family::FunctionBasedProblemFamily, ndim::Int)

Construct search space for `FunctionBasedProblem` with the given number of dimensions.
"""
function instantiate_search_space(family::FunctionBasedProblemFamily, ndim::Int)
    ndim >= numdims(family.reserved_ss) ||
        throw(ArgumentError("Cannot create $ndim-problem: number of dimensions less than reserved dimensions"))
    vcat(family.reserved_ss, symmetric_search_space(ndim - numdims(family.reserved_ss), family.range_per_dim))
end

"""
    instantiate(family::FunctionBasedProblemFamily, ndim::Int)

Construct `FunctionBasedProblem` with the given number of dimensions.
"""
instantiate(family::FunctionBasedProblemFamily, ndim::Int) =
    FunctionBasedProblem(family.objfunc, family.name, family.fitness_scheme,
                         instantiate_search_space(family, ndim), family.opt_value)

function instantiate(prob::OptimizationProblem, ndim::Int)
    numdims(prob) == ndim ||
        throw("Problem dimensions ($(numdims(prob))) and the requested dimension ($ndim) do not match")
    return prob # as is
end

@deprecate fixed_dim_problem instantiate

MinimizationProblemFamily(f::Function, name::String, range::ParamBounds, fmin::Float64) =
    FunctionBasedProblemFamily(f, name, MinimizingFitnessScheme, range, fmin)

MinimizationProblemFamily(f::Function, name::String, range::ParamBounds) =
    FunctionBasedProblemFamily(f, name, MinimizingFitnessScheme, range)

minimization_problem(f::Function, name::String, range::ParamBounds, ndim::Int) =
    instantiate(MinimizationProblemFamily(f, name, range), ndim)

minimization_problem(f::Function, range::ParamBounds, ndim::Int) =
    instantiate(MinimizationProblemFamily(f, "<unknown>", range), ndim)

minimization_problem(f::Function, name::String, range::ParamBounds, ndim::Int, fmin::Float64) =
    instantiate(MinimizationProblemFamily(f, name, range, fmin), ndim)

"""
    problem_set(ps::Dict{Any, FunctionBasedProblemFamily}, dims::Union{Int,Vector{Int}})

Construct a fixed-dimensional version of each problem from `ps`
for each dimension given in `dims`.

Returns a dictionary of problems.
"""
function problem_set(ps::Dict{Any, FunctionBasedProblemFamily}, dims::Vector{Int})
    next_free_index = 1
    # FIXME would using vector be more straightforward?
    result = Dict{Int, FunctionBasedProblem}()
    for d in dims
        for p in values(ps)
            result[next_free_index] = instantiate(p, d)
            next_free_index += 1
        end
    end
    return result
end

problem_set(ps::Dict{Any, FunctionBasedProblemFamily}, dim::Int) = problem_set(ps, [dim])
