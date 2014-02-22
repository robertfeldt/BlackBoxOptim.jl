# We have different sources for problem functions:
#  S1 = CEC 2013 competition on large-scale optimization
#  S2 = JADE paper http://150.214.190.154/EAMHCO/pdf/JADE.pdf
#  S3 = "Test Suite for the Special Issue of Soft Computing on Scalability of 
#        Evolutionary Algorithms and other Metaheuristics for Large Scale 
#        Continuous Optimization Problems", http://sci2s.ugr.es/eamhco/functions1-19.pdf
# Our primary focus is to implement all the problems from S1 since our
# focus is on large-scale optimization but these problems also can be used
# in lower dimensions.

#####################################################################
# Base functions.
#####################################################################
function sphere(x)
  sum(x.^2)
end

# Schwefel's ellipsoid.
function ellipsoid(x) 
  res = 0
  for(i in 1:length(x))
    res += sum(x[1:i])^2
  end
  res
end

function elliptic(x)
  D = length(x)
  condition = 1e+6
  coefficients = condition .^ linspace(0, 1, D)
  sum(coefficients .* x.^2)
end

function rastrigin(x)
  D = length(x)
  10 * D + sum( x.^2 ) - 10 * sum( cos( 2 * π * x ) )
end

function ackley(x)
  D = length(x)
  try
    20 - 20.*exp(-0.2.*sqrt(sum(x.^2)/D)) - exp(sum(cos(2 * π * x))/D) + e
  catch e
    # Sometimes we have gotten a DomainError from the cos function so we protect this call
    println(e)
    println("For input x = ", x)
    # Return a very large fitness value to indicate that this is NOT the solution we want.
    # TODO: Fix better way to handle this!
    9.99e100
  end
end

function schwefel1_2(x)
  D = length(x)
  partsums = zeros(D)
  partsum = 0
  for i in 1:D
    partsum += x[i]
    partsums[i] = partsum
  end
  sum(partsums.^2)
end

function rosenbrock(x)
  n = length(x)
  return( sum( 100*( x[2:n] - x[1:(n-1)].^2 ).^2 + ( x[1:(n-1)] - 1 ).^2 ) )
end

function step(x)
  sum(ceil(x + 0.5))
end

function griewank(x)
  n = length(x)
  1 + (1/4000)*sum(x.^2) - prod(cos(x ./ sqrt(1:n)))
end

function schwefel2_22(x)
  ax = abs(x)
  sum(ax) + prod(ax)
end

function schwefel2_21(x)
  maximum(abs(x))
end

# I'm unsure about this one since it does not return the expected minima at
# [1.0, 1.0].
function schwefel2_26(x)
  D = length(x)
  418.98288727243369 * D - sum(x .* sin(sqrt(abs(x))))
end

function cigar(x)
  x[1]^2 + 1e6 * sum(x[2:end].^2)
end

function cigtab(x)
  x[1]^2 + 1e8 * x[end]^2 + 1e4 * sum(x[2:(end-1)].^2)
end

# Shekel10 is a 4D, multi-minima, non-separable test problem. Our implementation
# is based on the C code in:
#   http://www.math.ntu.edu.tw/~wwang/cola_lab/test_problems/multiple_opt/multiopt_prob/Shekel10/Shekel10.c
Shekel10_A = [4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7; 2 9 2 9; 5 5 3 3; 8 1 8 1; 6 2 6 2; 7 3.6 7 3.6]
Shekel10_C = [0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5]
function shekel(x, a, c)
  sum = 0.0
  for i in 1:length(c)
    den = 0.0
    for j in 1:size(a, 2)
      den += (x[j] - a[i,j])^2
    end
    sum = sum - 1 / (den + c[i])
  end
  return sum
end
shekel10(x) = shekel(x, Shekel10_A, Shekel10_C)
Shekel10 = fixeddim_problem(shekel10; range = (0.0, 10.0), dims = 4, name = "Shekel10", fmins = [-10.5364])

# Shekel7 is a 4D, multi-minima, non-separable test problem. Our implementation
# is based on the C code in:
#   http://www.math.ntu.edu.tw/~wwang/cola_lab/test_problems/multiple_opt/multiopt_prob/Shekel7/Shekel7.c
Shekel7_A = [4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7; 2 9 2 9; 5 5 3 3]
Shekel7_C = [0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3]
shekel7(x) = shekel(x, Shekel7_A, Shekel7_C)
Shekel7 = fixeddim_problem(shekel7; range = (0.0, 10.0), dims = 4, name = "Shekel7", fmins = [-10.4029])

# Shekel5 is a 4D, multi-minima, non-separable test problem. Our implementation
# is based on the C code in:
#   http://www.math.ntu.edu.tw/~wwang/cola_lab/test_problems/multiple_opt/multiopt_prob/Shekel5/Shekel5.c
Shekel5_A = [4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7]
Shekel5_C = [0.1, 0.2, 0.2, 0.4, 0.4]
shekel5(x) = shekel(x, Shekel5_A, Shekel5_C)
Shekel5 = fixeddim_problem(shekel5; range = (0.0, 10.0), dims = 4, name = "Shekel5", fmins = [-10.1532])

# Hartman 6D is a multi-minima, non-separable test problem. Our implementation is based on:
#  http://www.sfu.ca/~ssurjano/hart6.html
Hartman6_alpha = [1.0 1.2 3.0 3.2]
Hartman6_A = [10 3 17 3.50 1.7 8; 0.05 10 17 0.1 8 14; 3 3.5 1.7 10 17 8; 17 8 0.05 10 0.1 14]
Hartman6_P = 1e-4 * [1312 1696 5569 124 8283 5886; 2329 4135 8307 3736 1004 9991; 2348 1451 3522 2883 3047 6650; 4047 8828 8732 5743 1091 381]
function hartman(x, alpha, A, P)
  sum = 0.0
  for i in 1:length(alpha)
    isum = 0.0
    for j in 1:size(A, 2)
      isum += A[i,j] * (x[j] - P[i,j])^2
    end
    sum -= alpha[i] * exp(-isum)
  end
  sum
end
hartman6(x) = hartman(x, Hartman6_alpha, Hartman6_A, Hartman6_P)
Hartman6 = fixeddim_problem(hartman6; range = (0.0, 1.0), dims = 6, name = "Hartman6", fmins = [-3.32237])

# Hartman 3D is a multi-minima, non-separable test problem. Our implementation is based on:
#  http://www.sfu.ca/~ssurjano/hart3.html
# However, we get a different global minima than the one stated on that page.
Hartman3_alpha = [1.0 1.2 3.0 3.2]
Hartman3_A = [3.0 10 30; 0.1 10 35; 3.0 10 30; 0.1 10 36]
Hartman3_P = 1e-4 * [3689 1170 2673; 4699 4387 7470; 1091 8732 5547; 381 5743 8828]
hartman3(x) = hartman(x, Hartman3_alpha, Hartman3_A, Hartman3_P)
# The minima should be -3.86278 but we get a different one so use that in the problem spec:
Hartman3 = fixeddim_problem(hartman3; range = (0.0, 1.0), dims = 3, name = "Hartman3", fmins = [-3.860038442])

#####################################################################
# S2 functions in addition to the base functions above. As stated
# in table II of the JADE paper: http://150.214.190.154/EAMHCO/pdf/JADE.pdf
#####################################################################

function quartic(x)
  D = length(x)
  sum( (1:D) .* x.^4 )
end

function noisy_quartic(x)
  quartic(x) + rand()
end

function s2_step(x)
  sum( ceil(x + 0.5).^2 )
end

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


#####################################################################
# Misc other interesting optimization functions and families.
#####################################################################

# This is a generator for the family of deceptive functions from the 
# Cuccu2011 paper on novelty-based restarts. We have vectorized it to allow
# more than 1D versions. The Cuccu2011 paper uses the following values for
# (l, w) = [(5, 0),  (15, 0),  (30, 0), 
#           (5, 2),  (15, 2),  (30, 2), 
#           (5, 10), (15, 10), (30, 10)]
# and notes that (15, 2) and (30, 2) are the most difficult instances.
function deceptive_cuccu2011(l, w)
  (x) -> begin
    absx = abs(x)
    sumabsx = sum(absx)
    if sumabsx <= 1
      return sum(x.^2)
    elseif sumabsx >= l+1
      return sum((absx - l).^2)
    else
      return (1 - 0.5 * sum(sin( (π * w * (absx - 1)) / l ).^2))
    end
  end
end

# Deceptive/hardest instances:
deceptive_cuccu2011_15_2 = deceptive_cuccu2011(15, 2)
deceptive_cuccu2011_30_2 = deceptive_cuccu2011(30, 2)

# Schwefel 2.13 is a hard one...
#function f=schwefel_213(x)
#global initial_flag
#persistent a b A alpha
#[ps,D]=size(x);
#if initial_flag==0
#    initial_flag=1;
#    load schwefel_213_data
#    if length(alpha)>=D
#        alpha=alpha(1:D);a=a(1:D,1:D);b=b(1:D,1:D);
#    else
#        alpha=-3+6*rand(1,D);
#        a=round(-100+200.*rand(D,D));
#        b=round(-100+200.*rand(D,D));
#    end
#    alpha=repmat(alpha,D,1);
#    A=sum(a.*sin(alpha)+b.*cos(alpha),2);
#end
#for i=1:ps
#    xx=repmat(x(i,:),D,1);
#    B=sum(a.*sin(xx)+b.*cos(xx),2);
#    f(i,1)=sum((A-B).^2,1);
#end

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
