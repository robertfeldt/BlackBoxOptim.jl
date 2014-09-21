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
Shekel10 = fixeddim_problem(shekel10; range = (0.0, 10.0), dims = 4, name = "Shekel10", fmins = [-10.5364])
Shekel7 = fixeddim_problem(shekel7; range = (0.0, 10.0), dims = 4, name = "Shekel7", fmins = [-10.4029])
Shekel5 = fixeddim_problem(shekel5; range = (0.0, 10.0), dims = 4, name = "Shekel5", fmins = [-10.1532])
Hartman6 = fixeddim_problem(hartman6; range = (0.0, 1.0), dims = 6, name = "Hartman6", fmins = [-3.32237])
Hartman3 = fixeddim_problem(hartman3; range = (0.0, 1.0), dims = 3, name = "Hartman3", fmins = [-3.860038442])

# We skip (for now) f12 and f13 in the JADE paper since they are penalized 
# functions which are quite nonstandard. We also skip f8 since we are unsure
# about its proper implementation.
JadeFunctionSet = {
  1   => anydim_problem("Sphere",        sphere,        (-100.0, 100.0), 0.0),
  2   => anydim_problem("Schwefel2.22",  schwefel2_22,  ( -10.0,  10.0), 0.0),
  3   => anydim_problem("Schwefel1.2",   schwefel1_2,   (-100.0, 100.0), 0.0),
  4   => anydim_problem("Schwefel2.21",  schwefel2_21,  (-100.0, 100.0), 0.0),
  5   => anydim_problem("Rosenbrock",    rosenbrock,    ( -30.0,  30.0), 0.0),
  6   => anydim_problem("Step",          s2_step,       (-100.0, 100.0), 0.0),
  7   => anydim_problem("Noisy quartic", noisy_quartic, ( -30.0,  30.0)),
#  8   => anydim_problem("Schwefel2.26",  schwefel2_26,  (-500.0, 500.0)),
  9   => anydim_problem("Rastrigin",     rastrigin,     ( -5.12,  5.12), 0.0),
  10  => anydim_problem("Ackley",        ackley,        ( -32.0,  32.0), 0.0),
  11  => anydim_problem("Griewank",      griewank,      (-600.0, 600.0), 0.0)
}


#####################################################################
# S3 Base functions.
#####################################################################


#####################################################################
# S3 Transformations
#####################################################################

# A TransformedProblem just makes a few changes in a sub-problem but refers
# most func calls to it. Concrete types must implement a sub_problem func.
abstract TransformedProblem <: OptimizationProblem
search_space(tp::TransformedProblem) = search_space(sub_problem(tp))
is_fixed_dimensional(tp::TransformedProblem) = is_fixed_dimensional(sub_problem(tp))
numfuncs(tp::TransformedProblem) = numfuncs(sub_problem(tp))
numdims(tp::TransformedProblem) = numdims(sub_problem(tp))
fmins(tp::TransformedProblem) = fmins(sub_problem(tp))
name(tp::TransformedProblem) = name(sub_problem(tp))

# A ShiftedAndBiasedProblem shifts the minimum value and biases the returned 
# function values.
type ShiftedAndBiasedProblem <: TransformedProblem
  xshift::Array{Float64, 1}
  funcshift::Float64
  subp::OptimizationProblem

  ShiftedAndBiasedProblem(sub_problem::OptimizationProblem; 
    xshift = false, funcshift = 0.0) = begin
    xshift = (xshift != false) ? xshift : rand_individual(search_space(sub_problem))
    new(xshift[:], funcshift, sub_problem)
  end
end

sub_problem(sp::ShiftedAndBiasedProblem) = sp.subp

is_fixed_dimensional(p::ShiftedAndBiasedProblem) = is_fixed_dimensional(sub_problem(p))

# Evaluate by first shifting x and then biasing the returned function value.
evalfunc(x, i::Int64, sp::ShiftedAndBiasedProblem) = begin
  ofunc(sub_problem(sp), i)(x - sp.xshift) + sp.funcshift
end

shifted(p::OptimizationProblem; funcshift = 0.0) = ShiftedAndBiasedProblem(p; 
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

example_problems = {
  "Sphere" => JadeFunctionSet[1],
  "Rosenbrock" => JadeFunctionSet[5],
  "Schwefel2.22" => JadeFunctionSet[2],
  "Schwefel1.2" => JadeFunctionSet[3],
  "Schwefel2.21" => JadeFunctionSet[4],
  "Step" => JadeFunctionSet[6],
  "Rastrigin" => JadeFunctionSet[9],
  "Ackley" => JadeFunctionSet[10],
  "Griewank" => JadeFunctionSet[11],
  "Ellipsoid" => anydim_problem("Ellipsoid", ellipsoid, (-65.536, 65.536), 0.0),
  "Cigar" => anydim_problem("Cigar", cigar, (-100.0, 100.0), 0.0),
  "DeceptiveCuccu2011_15_2" => anydim_problem("DeceptiveCuccu2011_15_2", deceptive_cuccu2011_15_2, (-100.0, 100.0), 0.0),
  "Shekel10" => Shekel10,
  "Shekel7" => Shekel7,
  "Shekel5" => Shekel5,
  "Hartman6" => Hartman6,
  "Hartman3" => Hartman3
}
