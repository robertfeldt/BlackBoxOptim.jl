using NLopt
using DataFrames

function rosenbrock(x)
  n = length(x)
  return( sum( 100*( x[2:n] - x[1:(n-1)].^2 ).^2 + ( x[1:(n-1)] - 1 ).^2 ) )
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

numfevals = 0 # keep track of # function evaluations

function make_objfunc(funcToOptimize)
  function objf(x::Vector, grad::Vector)
    global numfevals
    numfevals::Int += 1
    funcToOptimize(x)
  end
  return objf
end

function compare_nlopt_algs(funcToOpt, dim;
  LowBound = -10.0,
  HighBound = 10.0,
  NumReps = 3,
  MaxTimes = [1.0, 3.0, 10.0], 
  Algs = [:LN_COBYLA, :LN_BOBYQA, :LN_NEWUOA, :LN_PRAXIS, :LN_NELDERMEAD, :LN_SBPLX,
    :GN_DIRECT, :GN_DIRECT_L, :GN_CRS2_LM, :GN_ISRES, :GN_ESCH]
  )

  function objfunc(x::Vector, grad::Vector)
    # We don't have gradient...

    global numfevals
    numfevals::Int += 1

    funcToOpt(x)
  end

  res = DataFrame()

  for maxt in MaxTimes
    for alg in Algs
      for rep in 1:NumReps

        opt = Opt(alg, dim)
        lower_bounds!(opt, (LowBound  * ones(dim)))
        upper_bounds!(opt, (HighBound * ones(dim)))
        maxtime!(opt, maxt)
        min_objective!(opt, objfunc)

        println("Solving...")

        global numfevals = 0
        tic()
        (fmin, solution, ret) = optimize(opt, randn(dim))
        t = toq()

        println("Dims = ", dim)
        println("Alg = ", alg)
        println("time given = ", maxt, " seconds")
        println("time taken = ", t, " seconds (", round(t/maxt*100.0, 2), "%)")
        println("fmin = ", fmin)
        println("Fevals = ", numfevals)
        println("Return value = ", ret, "\n\n")

        res = rbind(res, DataFrame(Alg = alg, D = dim, Rep = rep, Fevals = numfevals,
          TimeGiven = maxt, TimeTaken = t, TimeUsed = round(100.0*t/maxt, 2), fmin = fmin))
      end
    end
  end
  res
end

r1 = compare_nlopt_algs(rosenbrock, 8; NumReps = 1, MaxTimes = [0.1, 1.0, 3.0])

numfevals = 0

function rosenbrock_objfunc(x::Vector, grad::Vector)
    # We don't have gradient...

  global numfevals
  numfevals::Int += 1

  rosenbrock(x)
end

# Given an array of NLopt algs and a vector of relative times to apply them
# we run them in sequence on the value of the previous one.
function combine_nlopt_algs(objfunc, dim;
  Algs = nothing, RelTimes = nothing,
  LowBound = -10.0,
  HighBound = 10.0,
  MaxTime = nothing
  )

  if Algs == nothing
    numpairs = iceil(log10(dim))
    Algs = vcat([:GN_DIRECT, :LN_PRAXIS], random_pairs_of_algs(numpairs), shuffle(LocalNLoptAlgs))
    RelTimes = vcat(rand(2 + 2*numpairs), 0.05 * rand(length(LocalNLoptAlgs)))
  end

  if RelTimes == nothing
    rand(length(Algs))
  end

  if MaxTime == nothing
    MaxTime = 3*sqrt(dim)
  end

  numsteps = length(Algs)
  @assert numsteps == length(RelTimes)

  # Ensure percentages
  percenttimes = RelTimes / sum(RelTimes)

  fmin = ret = 0
  global numfevals = 0

  xbest = x = randn(dim)
  fbest = objfunc(xbest, zeros(dim))

  for i in 1:numsteps
    opt = Opt(Algs[i], dim)
    lower_bounds!(opt, (LowBound  * ones(dim)))
    upper_bounds!(opt, (HighBound * ones(dim)))
    maxt = MaxTime * percenttimes[i]
    maxtime!(opt, maxt)
    ftol_rel!(opt, 1e-7)
    xtol_rel!(opt, 1e-7)
    min_objective!(opt, objfunc)
    (fmin, x, ret) = optimize(opt, x)
    println(Algs[i], " for ", round(maxt, 3), " s, fmin = ", fmin, ", ret = ", ret)
    if fmin < fbest
      fbest = fmin
      xbest = x
    else
      x = xbest
    end
  end

  (fbest, xbest, ret, copy(numfevals))
end

LocalNLoptAlgs = [:LN_COBYLA, :LN_BOBYQA, :LN_NEWUOA, :LN_PRAXIS, :LN_NELDERMEAD, :LN_SBPLX]
GlobalNLoptAlgs = [:GN_DIRECT, :GN_DIRECT_L, :GN_CRS2_LM, :GN_ISRES, :GN_ESCH]

# Randomly select np pairs of algs where the first is global and the second is local.
function random_pairs_of_algs(np)
  algs = Symbol[]
  for i in 1:np
    push!(algs, shuffle(GlobalNLoptAlgs)[1])
    push!(algs, shuffle(LocalNLoptAlgs)[1])
  end
  algs
end

numpairs = 4 
rosenbrock_of = make_objfunc(rosenbrock)
r = combine_nlopt_algs(rosenbrock_of, 64; MaxTime = 30.0, 
  Algs = vcat([:LN_PRAXIS], random_pairs_of_algs(numpairs)), RelTimes = rand(1+2*numpairs))

rastrigin_of = make_objfunc(rastrigin)
numpairs = 4 
r = combine_nlopt_algs(rastrigin_of, 8; MaxTime = 3.0, 
  Algs = vcat([:LN_PRAXIS], random_pairs_of_algs(numpairs)), RelTimes = rand(1+2*numpairs))

ackley_of = make_objfunc(ackley)
numpairs = 4 
r = combine_nlopt_algs(ackley_of, 8; MaxTime = 3.0, 
  Algs = vcat([:GN_DIRECT, :LN_PRAXIS], random_pairs_of_algs(numpairs)), RelTimes = rand(2+2*numpairs))
r = combine_nlopt_algs(rastrigin_of, 8; MaxTime = 3.0, 
  Algs = vcat([:GN_DIRECT, :LN_PRAXIS], random_pairs_of_algs(numpairs)), RelTimes = rand(2+2*numpairs))
r = combine_nlopt_algs(rosenbrock_of, 16; MaxTime = 3.0, 
  Algs = vcat([:GN_DIRECT, :LN_PRAXIS], random_pairs_of_algs(numpairs)), RelTimes = rand(2+2*numpairs))

r = combine_nlopt_algs(rosenbrock_of, 16)
r = combine_nlopt_algs(rosenbrock_of, 32)
r = combine_nlopt_algs(rosenbrock_of, 64)
r = combine_nlopt_algs(rosenbrock_of, 512)
r = combine_nlopt_algs(rastrigin_of, 512)
r = combine_nlopt_algs(ackley_of, 512)