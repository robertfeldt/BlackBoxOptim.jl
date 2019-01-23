# We have different sources for problem functions:
#  S1 = CEC 2013 competition on large-scale optimization
#  S2 = JADE paper http://150.214.190.154/EAMHCO/pdf/JADE.pdf
#  S3 = "Test Suite for the Special Issue of Soft Computing on Scalability of
#        Evolutionary Algorithms and other Metaheuristics for Large Scale
#        Continuous Optimization Problems", http://sci2s.ugr.es/eamhco/functions1-19.pdf
# Our primary focus is to implement all the problems from S1 since our
# focus is on large-scale optimization but these problems also can be used
# in lower dimensions.
include(joinpath(dirname(@__FILE__()), "single_objective_base_functions.jl"))
Shekel10 = minimization_problem(shekel10, "Shekel10", (0.0, 10.0), 4, -10.5364)
Shekel7 = minimization_problem(shekel7, "Shekel7", (0.0, 10.0), 4, -10.4029)
Shekel5 = minimization_problem(shekel5, "Shekel5", (0.0, 10.0), 4, -10.1532)
Hartman6 = minimization_problem(hartman6, "Hartman6", (0.0, 1.0), 6, -3.32237)
Hartman3 = minimization_problem(hartman3, "Hartman3", (0.0, 1.0), 3, -3.860038442)

"""
JADE collection of optimization problems.

We skip (for now) `f12` and `f13` in the JADE paper since they are penalized
functions which are quite nonstandard. We also skip `f8` since we are unsure
about its proper implementation.
"""
const JadeFunctionSet = Dict{Int,FunctionBasedProblemFamily}(
    1   => MinimizationProblemFamily(sphere, "Sphere", (-100.0, 100.0), 0.0),
    2   => MinimizationProblemFamily(schwefel2_22,  "Schwefel2.22",  ( -10.0,  10.0), 0.0),
    3   => MinimizationProblemFamily(schwefel1_2,   "Schwefel1.2",   (-100.0, 100.0), 0.0),
    4   => MinimizationProblemFamily(schwefel2_21,  "Schwefel2.21",  (-100.0, 100.0), 0.0),
    5   => MinimizationProblemFamily(rosenbrock,    "Rosenbrock",    ( -30.0,  30.0), 0.0),
    6   => MinimizationProblemFamily(s2_step,       "Step",          (-100.0, 100.0), 0.0),
    7   => MinimizationProblemFamily(noisy_quartic, "Noisy quartic", ( -30.0,  30.0)),
#    8   => MinimizationProblemFamily(schwefel2_26,  "Schwefel2.26",  (-500.0, 500.0)),
    9   => MinimizationProblemFamily(rastrigin,     "Rastrigin",     ( -5.12,  5.12), 0.0),
    10  => MinimizationProblemFamily(ackley,        "Ackley",        ( -32.0,  32.0), 0.0),
    11  => MinimizationProblemFamily(griewank,      "Griewank",      (-600.0, 600.0), 0.0)
)

#####################################################################
# S3 Base functions.
#####################################################################


#####################################################################
# S3 Transformations
#####################################################################

"""
`TransformedProblem{FS, P}` is an abstract class for optimization problems derived
from some original problem of type `P` by introducing a few changes.

The concrete derived types must implement:
* `objfunc()` method
* `fitness_scheme()` and `search_space()`
"""
abstract type TransformedProblem{FS, P <: OptimizationProblem} <: OptimizationProblem{FS} end

orig_problem(tp::TransformedProblem) = tp.orig_problem # default implementation
orig_problem_type(tp::Type{<:TransformedProblem{P}}) where P = P
name(tp::TransformedProblem) = "TransformedProblem("*name(sub_problem(tp))*")"

"""
A `TransformedProblem` subclass that shifts the objective argument and offsets the value:
``f(x - xshift) + fitshift`` → min/max, ``x`` ∈ ``X`` (``X`` stays intact).
"""
struct ShiftedAndBiasedProblem{FS, F, P <: OptimizationProblem} <: TransformedProblem{FS, P}
    orig_problem::P

    xshift::Vector{Float64}
    fitshift::F
    opt_value::F

    function ShiftedAndBiasedProblem(orig_problem::P;
        xshift = nothing,
        fitshift = nothing,
        opt_value = nothing
    ) where P <: OptimizationProblem{FS} where FS
        (xshift === nothing) && (xshift = fill(0.0, numdims(orig_problem)))
        F = fitness_type(FS)
        (fitshift === nothing) && (fitshift = zero(F)) # FIXME zerofitness() to support tuples
        # derive the default optimal value,
        # which would be not correct if `xshift` moves the optimum in/out of the search space
        (opt_value === nothing) && (opt_value = BlackBoxOptim.opt_value(orig_problem) + fitshift)
        new{FS, F, P}(orig_problem, copy(xshift), fitshift, opt_value)
    end
end

name(p::ShiftedAndBiasedProblem) = "ShiftedAndBiased("*name(orig_problem(p))*")"
fitness_scheme(tp::TransformedProblem) = fitness_scheme(orig_problem(tp))
search_space(tp::TransformedProblem) = search_space(orig_problem(tp)) # the space stays the same
opt_value(tp::TransformedProblem) = tp.opt_value

fitness(x, sp::ShiftedAndBiasedProblem) =
    fitness(x - sp.xshift, orig_problem(sp)) + sp.fitshift

#####################################################################
# S1 Base functions. Typically slightly transformed to break symmetry
#   and introduce irregularities.
#####################################################################
s1_sphere = sphere

function s1_elliptic(x)
    xt = t_irreg(x)
    elliptic(xt)
end

function s1_rastrigin(x)
    xt = t_diag(t_asy(t_irreg(x), 0.2), 10)
    rastrigin(xt)
end

function s1_ackley(x)
    xt = t_diag(t_asy(t_irreg(x), 0.2), 10)
    ackley(xt)
end

function s1_schwefel(x)
    xt = t_asy(t_irreg(x), 0.2)
    schwefel1_2(xt)
end

s1_rosenbrock = rosenbrock


#####################################################################
# S1 Transformations
#####################################################################

"""
Transform symmetric `f` into asymmetric objective function.
"""
function t_asy(f, beta)
    D = length(f)
    g = copy(f)
    temp = beta * range(0, stop=1, length=D)
    ind = collect(1:D)[f .> 0]
    t = f[ind] .^ (1 + temp[ind] .* sqrt(f[ind]))
    setindex!(g, t, ind)
    return g
end

"""
Transform `f` into objective function with ill-conditioning effect.
"""
function t_diag(f, alpha)
    D = length(f)
    scales = sqrt(alpha) .^ range(0, stop=1, length=D)
    return scales .* f
end

"""
Transform `f` into objective function with smooth local irregularities.
"""
function t_irreg(f)
    a = 0.1
    g = copy(f)
    indices = collect(1:length(f))

    idxp = indices[f .> 0]
    t = log(f[idxp])/a
    r = exp(t + 0.49*(sin(t) + sin(0.79*t))).^a
    setindex!(g, r, idxp)

    idxn = indices[f .< 0]
    t = log(-f[idxn])/a
    r = -exp(t + 0.49*(sin(0.55*t) + sin(0.31*t))).^a
    setindex!(g, r, idxn)

    return g
end

function xshifted(n, f)
    move = 10.0 * randn(n, 1)
    transformed_f(x) = f(x .- move)
end

function xrotatedandshifted(n, f, shiftAmplitude = 1.0, rotateAmplitude = 1.0)
    shift = shiftAmplitude * randn(n, 1)
    rotmatrix = rotateAmplitude * rand(n, n)
    transformed_f(x) = f(rotmatrix * (x .- shift))
end

const example_problems = Dict{String, Union{OptimizationProblem,FunctionBasedProblemFamily}}(
    "Sphere" => JadeFunctionSet[1],
    "Rosenbrock" => JadeFunctionSet[5],
    "Schwefel2.22" => JadeFunctionSet[2],
    "Schwefel1.2" => JadeFunctionSet[3],
    "Schwefel2.21" => JadeFunctionSet[4],
    "Step" => JadeFunctionSet[6],
    "Rastrigin" => JadeFunctionSet[9],
    "Ackley" => JadeFunctionSet[10],
    "Griewank" => JadeFunctionSet[11],
    "Ellipsoid" => MinimizationProblemFamily(ellipsoid, "Ellipsoid", (-65.536, 65.536), 0.0),
    "Cigar" => MinimizationProblemFamily(cigar, "Cigar", (-100.0, 100.0), 0.0),
    "DeceptiveCuccu2011_15_2" => MinimizationProblemFamily(deceptive_cuccu2011_15_2, "DeceptiveCuccu2011_15_2", (-100.0, 100.0), 0.0),
    "Shekel10" => Shekel10,
    "Shekel7" => Shekel7,
    "Shekel5" => Shekel5,
    "Hartman6" => Hartman6,
    "Hartman3" => Hartman3,
    "Tsallis1996" => MinimizationProblemFamily(energy_tsallis1996, "Tsallis1996", (-100.0, 100.0), 0.0),
)
