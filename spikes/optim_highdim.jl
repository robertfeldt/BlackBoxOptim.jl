using Optim
using DataFrames

include(joinpath("..", "src", "problems", "single_objective_base_functions.jl"))

function xrotatedandshifted(n, f, shiftAmplitude = 1.0, rotateAmplitude = 1.0)
  shift = shiftAmplitude * randn(n, 1)
  rotmatrix = rotateAmplitude * rand(n, n)
  transformed_f(x) = f(rotmatrix * (x .- shift))
end

OptimAlgs = [:bfgs, :nelder_mead, :l_bfgs, :cg,
  :gradient_descent, :momentum_gradient_descent, :simulated_annealing]

function compare_algs(funcToOpt, dim;
  LowBound = -10.0,
  HighBound = 10.0,
  NumReps = 3,
  NumIterations = 2000,
  TraceMode = :compact,
  Algs = OptimAlgs)

  global tf

  function objfunc(x::Vector)
    global numfevals
    numfevals::Int += 1
    #funcToOpt(x)
    tf(x)
  end

  dfs = DataFrame[]

  for alg in Algs
    for rep in 1:NumReps
      global tf = xrotatedandshifted(dim, funcToOpt)
      global numfevals = 0
      if TraceMode != :silent
        println("\nSolving...")
      end
      try
        tic()
        res = optimize(objfunc, randn(dim), method = alg, iterations = NumIterations)
        t = toq()
        if TraceMode != :silent
          println("  Dims = ", dim)
          println("  Alg = ", alg)
          println("  time taken = ", t, " seconds")
          println("  fmin = ", res.f_minimum)
          println("  Fevals = ", numfevals)
        end
        push!(dfs, DataFrame(Alg = alg, D = dim, Rep = rep, Iterations = NumIterations,
          Fevals = numfevals, TimeTaken = t, fmin = res.f_minimum))
      catch err
        println("ERROR: ", err)
      end
    end
  end
  vcat(dfs...)
end

# Short one so things are compiled
compare_algs(rastrigin, 2; NumReps = 1, TraceMode = :silent);

# They use a very uneven num of fevals though...
BestAlgs = [:bfgs, :nelder_mead, :simulated_annealing]

r1 = compare_algs(rastrigin, 100; NumReps = 3, NumIterations = 100000,
  Algs = BestAlgs)
println(r1)
x = by(r1, [:Alg, :D], df -> DataFrame(
  N = size(df, 1),
  mFevals = median(df[:,:Fevals]),
  mTime = median(df[:,:TimeTaken]),
  mFmin = median(df[:,:fmin])))
println(x)
sort!(x, cols = [:mFmin])
println(x)
