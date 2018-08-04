using NLopt
using DataFrames

include(joinpath("..", "src", "problems", "single_objective_base_functions.jl"))

abstract type Evaluator end

mutable struct NonGradientEvaluator <: Evaluator
  func::Function
  fevals::Int
  NonGradientEvaluator(func) = new(func, 0)
end
import Base.reset
reset(e::Evaluator) = (e.fevals = 0)
fevals(e::Evaluator) = e.fevals

function makeobjfunc(e::Evaluator)
  ofunc = (x::Vector, grad::Vector) -> begin
    e.fevals::Int += 1
    e.func(x)
  end
  ofunc
end

LocalNLoptAlgs = [:LN_COBYLA, :LN_BOBYQA, :LN_NEWUOA, :LN_PRAXIS, :LN_NELDERMEAD, :LN_SBPLX]
GlobalNLoptAlgs = [:GN_DIRECT, :GN_DIRECT_L, :GN_ORIG_DIRECT, :GN_ORIG_DIRECT_L, 
  :GN_CRS2_LM, :GN_ISRES, :GN_ESCH]

# Randomly select np pairs of algs where the first is global and the second is local.
function random_pairs_of_algs(np)
  algs = Symbol[]
  for i in 1:np
    push!(algs, shuffle(GlobalNLoptAlgs)[1])
    push!(algs, shuffle(LocalNLoptAlgs)[1])
  end
  algs
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
        println("time taken = ", t, " seconds (", round(t/maxt*100.0, digits=2), "%)")
        println("fmin = ", fmin)
        println("Fevals = ", numfevals)
        println("Return value = ", ret, "\n\n")

        res = rbind(res, DataFrame(Alg = alg, D = dim, Rep = rep, Fevals = numfevals,
          TimeGiven = maxt, TimeTaken = t, TimeUsed = round(100.0*t/maxt, digits=2), fmin = fmin))
      end
    end
  end
  res
end

# Given an array of NLopt algs and a vector of relative times to apply them
# we run them in sequence on the value of the previous one.
function combine_nlopt_algs(func, dim;
  Algs = nothing, RelTimes = nothing,
  LowBound = -100.0,
  HighBound = 100.0,
  MaxTime = nothing,
  ResetToBestProbability = 0.80,
  FTOL_REL = -Inf,
  FTOL_ABS = -Inf
  )

  if Algs == nothing
    numpairs = iceil(log10(dim))
    Algs = vcat([:GN_DIRECT_L, :LN_PRAXIS], random_pairs_of_algs(numpairs), 
      shuffle(LocalNLoptAlgs), [:LN_PRAXIS])
    RelTimes = vcat([1, 3], rand(1:3, 2*numpairs), 0.30 * rand(length(LocalNLoptAlgs)), 2)
  end

  if RelTimes == nothing
    RelTimes = rand(length(Algs))
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

  evaluator = NonGradientEvaluator(func)
  objfunc = makeobjfunc(evaluator)

  xbest = x = randn(dim)
  fbest = objfunc(xbest, zeros(dim))
  numtolreached = fevalsbest = 0

  start_time = time()

  println("\nStarting solve sequence...")
  for i in 1:numsteps
    opt = Opt(Algs[i], dim)
    lower_bounds!(opt, (LowBound  * ones(dim)))
    upper_bounds!(opt, (HighBound * ones(dim)))
    maxt = MaxTime * percenttimes[i]
    maxtime!(opt, maxt)
    ftol_rel!(opt, FTOL_REL)
    ftol_abs!(opt, FTOL_ABS)
    min_objective!(opt, objfunc)
    try
      tic()
      (fmin, x, ret) = optimize(opt, x)
      t = toq()
      println(Algs[i], " for ", round(t, digits=3), " s (", round(t/maxt*100.0, digits=2), "%), fmin = ", fmin, ", ret = ", ret)
    catch err
      println(Algs[i], " failed. Skipping...")
    end
    nowtime = time()
    if fmin < fbest
      fbest = fmin
      xbest = x
      fevalsbest = fevals(evaluator)
      timebest = nowtime
    else
      #if rand() < ResetToBestProbability
      #  x = xbest
      #end
    end
    elapsed = nowtime - start_time
    if (nowtime - timebest) > 0.20 * MaxTime
      if in(ret, [:SUCCESS, :FTOL_REACHED, :XTOL_REACHED, :ROUNDOFF_LIMITED])
        if numtolreached < 1
          # If first time to converge we try once more from a new random position
          numtolreached += 1
          x = randn(dim)
          println("Restart since ret = ", ret)
        else
          println("Return since ret = ", ret, " and numtolreached = ", numtolreached)
          return (fbest, xbest, :CONVERGED, fevalsbest, fevalsbest/dim, fevals(evaluator), fevals(evaluator)/dim, time() - start_time)
        end
      end
    end
  end

  return (fbest, xbest, ret, fevalsbest, fevalsbest/dim, fevals(evaluator), fevals(evaluator)/dim, time() - start_time)
end

dfmapper(r) = begin
  DataFrame(nfevalsbest = r[4], nfevals = r[6], nfbpd = r[5], fmin = r[1])
end

SOPS = [sphere, ellipsoid, elliptic, cigar, cigtab, ackley, schwefel1_2, 
  schwefel2_22, schwefel2_21, rastrigin, rosenbrock]

SOPSNonzero = [step, griewank, schwefel2_26, shekel10]

function xrotatedandshifted(n, f, shiftAmplitude = 1.0, rotateAmplitude = 1.0)
  shift = shiftAmplitude * randn(n, 1)
  rotmatrix = rotateAmplitude * rand(n, n)
  transformed_f(x) = f(rotmatrix * (x .- shift))
end

#res8 = map((sop) -> combine_nlopt_algs(sop, 8), SOPS)
#resrot8 = map((sop) -> combine_nlopt_algs(xrotatedandshifted(8, sop), 8), SOPS)

#res64 = map((sop) -> combine_nlopt_algs(sop, 64), SOPS)
#resrot64 = map((sop) -> combine_nlopt_algs(xrotatedandshifted(64, sop), 64), SOPS)

#res256 = map((sop) -> combine_nlopt_algs(sop, 256), SOPS)
#resrot256 = map((sop) -> combine_nlopt_algs(xrotatedandshifted(256, sop), 256), SOPS)

function compare_optimizers(optimizers, problems, dim, numreps = 10;
  dataframe_mapper = identity)

  timestamp = strftime("%Y%m%d_%H%M%S", time())
  filename = "results_nlopt_$(timestamp).csv"
  res = DataFrame()
  its = 0
  for rep in 1:numreps
    for (probname, origproblem) in problems
      problem = xrotatedandshifted(dim, origproblem)
      for (optname, opt) in optimizers
        its += 1
        println("Optimization ", its)
        tic()
        r = opt(problem, dim)
        t = toq()
        res = rbind(res, cbind(
          DataFrame(Opt = optname, Prob = probname, Dim = dim, Time = round(t, digits=2)),
          dataframe_mapper(r)
          ))
        writetable(filename, res)
        println("Wrote results to file: ", filename)
      end
    end
  end
  return res, filename
end

Problems = {
  "Sphere      " => sphere,
  "Rosenbrock  " => rosenbrock,
  "Schwefel2_22" => schwefel2_22,
  "Schwefel2_21" => schwefel2_21, 
  "Rastrigin   " => rastrigin,
}

DP = [:GN_DIRECT_L, :LN_PRAXIS]
DPN = [:GN_DIRECT_L, :LN_PRAXIS, :LN_NELDERMEAD]

Optimizers = {
  "DirL1+Prax3       " => (pr, d) -> combine_nlopt_algs(pr, d; Algs = DP, RelTimes = [1,3]),
  "DirL1+Prax2+NelM2 " => (pr, d) -> combine_nlopt_algs(pr, d; Algs = DPN, RelTimes = [1,2,2]),
  "RandNLopt non-grad" => (pr, d) -> combine_nlopt_algs(pr, d),
}

function compare_and_save(opts, probs, dims, reps, mapper)
  for dim in dims
    res, filename = compare_optimizers(opts, probs, dim, reps; dataframe_mapper = mapper)
    println("Wrote results to file: ", filename)
  end
end

compare_and_save(Optimizers, Problems, [256, 1024], 1, dfmapper)