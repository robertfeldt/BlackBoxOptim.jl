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

# We skip (for now) f12 and f13 in the JADE paper since they are penalized
# functions which are quite nonstandard. We also skip f8 since we are unsure
# about its proper implementation.
const JadeFunctionSet = @compat Dict{Int,FunctionBasedProblemFamily}(
  1   => MinimizationProblemFamily(sphere, "Sphere", (-100.0, 100.0), 0.0),
  2   => MinimizationProblemFamily(schwefel2_22,  "Schwefel2.22",  ( -10.0,  10.0), 0.0),
  3   => MinimizationProblemFamily(schwefel1_2,   "Schwefel1.2",   (-100.0, 100.0), 0.0),
  4   => MinimizationProblemFamily(schwefel2_21,  "Schwefel2.21",  (-100.0, 100.0), 0.0),
  5   => MinimizationProblemFamily(rosenbrock,    "Rosenbrock",    ( -30.0,  30.0), 0.0),
  6   => MinimizationProblemFamily(s2_step,       "Step",          (-100.0, 100.0), 0.0),
  7   => MinimizationProblemFamily(noisy_quartic, "Noisy quartic", ( -30.0,  30.0)),
#  8   => MinimizationProblemFamily(schwefel2_26,  "Schwefel2.26",  (-500.0, 500.0)),
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

# A TransformedProblem just makes a few changes in a sub-problem but refers
# most func calls to it. Concrete types must implement a sub_problem func.
abstract TransformedProblem{FS<:FitnessScheme} <: OptimizationProblem{FS}
search_space(tp::TransformedProblem) = search_space(sub_problem(tp))
is_fixed_dimensional(tp::TransformedProblem) = is_fixed_dimensional(sub_problem(tp))
numfuncs(tp::TransformedProblem) = numfuncs(sub_problem(tp))
numdims(tp::TransformedProblem) = numdims(sub_problem(tp))
fmins(tp::TransformedProblem) = fmins(sub_problem(tp))
name(tp::TransformedProblem) = name(sub_problem(tp))

# A ShiftedAndBiasedProblem shifts the minimum value and biases the returned
# function values.
type ShiftedAndBiasedProblem{FS<:FitnessScheme} <: TransformedProblem{FS}
  xshift::Array{Float64, 1}
  funcshift::Float64
  subp::OptimizationProblem{FS}

  ShiftedAndBiasedProblem(sub_problem::OptimizationProblem{FS};
    xshift = false, funcshift = 0.0) = begin
    xshift = (xshift != false) ? xshift : rand_individual(search_space(sub_problem))
    new(xshift[:], funcshift, sub_problem)
  end
end

sub_problem(sp::ShiftedAndBiasedProblem) = sp.subp

is_fixed_dimensional(p::ShiftedAndBiasedProblem) = is_fixed_dimensional(sub_problem(p))

# Evaluate by first shifting x and then biasing the returned function value.
evalfunc(x, i, sp::ShiftedAndBiasedProblem) = begin
  ofunc(sub_problem(sp), i)(x - sp.xshift) + sp.funcshift
end

shifted{FS<:FitnessScheme}(p::OptimizationProblem{FS}; funcshift = 0.0) = ShiftedAndBiasedProblem{FS}(p;
  funcshift = funcshift)


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

# This transformation function is used to break the symmetry of symmetric
# functions.
function t_asy(f, beta)
  D = length(f)
  g = copy(f)
  temp = beta * linspace(0, 1, D)
  ind = collect(1:D)[f .> 0]
  t = f[ind] .^ (1 + temp[ind] .* sqrt(f[ind]))
  setindex!(g, t, ind)
  g
end

# This transformation is used to create the ill-conditioning effect.
function t_diag(f, alpha)
  D = length(f)
  scales = sqrt(alpha) .^ linspace(0, 1, D)
  scales .* f
end

# This transformation is used to create smooth local irregularities.
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

   g
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

const example_problems = @compat Dict{String,Any}( #FIXME use Union{Optimization,FunctionBasedProblemFamily}
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
