function sumstats(v, f = (x) -> @sprintf("%.2f", x))
  try
    v = convert(Array{Float64,1}, v)
  catch
  end
  "median = $(f(median(v))), mean = $(f(mean(v))) +/- $(f(std(v))), range = [$(f(minimum(v))), $(f(maximum(v)))]"
end

function format_time(t)
  if t < 5e-1
    @sprintf("%.2f ms", t*1e3)
  elseif t < 60.0
    @sprintf("%.2f s", t)
  elseif t < 30*60
    @sprintf("%.2f min", t/60.0)
  else
    @sprintf("%.2f hours", t/3600.0)
  end
end

function unique_filename(prefix = "result", suffix = ".txt")
  join([prefix, strftime("_%Y%m%d_%H%M%S_$(rand(1:int(1e6)))", time()), suffix])
end

function csvfile(header; fileprefix = "experiment", 
  filepath = unique_filename(fileprefix, ".csv"), sep = ",")
  fh = open(filepath, "w")
  println(fh, join(header, sep))
  fh
end

# Run repeated experimental runs of a search function for each problem
# in a problem set.
function repeated_runs(searchf, problem_set, num_runs = 10; experiment = "exp")
  problems = collect(problem_set)
  num_problems = length(problems)
  fevals = zeros(Int64, num_runs, num_problems)
  fbests = zeros(num_runs, num_problems)
  times = zeros(num_runs, num_problems)
  reason_counts = [Dict{ASCIIString, Int64}() for i in 1:num_problems]

  csvfh = csvfile(["Experiment", "Date", "Time", "RunId", 
    "Problem", "Dimensions", "TerminationReason", "FuncEvals", "ElapsedTime", "Fitness"]; fileprefix = experiment)

  for(i in 1:num_runs)
    # Random order of running each problem
    problem_indices = shuffle(collect(1:num_problems))
    for(pi in problem_indices)
      num, prob = problems[pi]
      optfunc = (x) -> BlackBoxOptim.eval1(x, prob)
      dims = BlackBoxOptim.numdims(prob)
      println("Run $(i) of problem $(prob.name)")
      start_time = time()
      tic()
      x, fbests[i,pi], fevals[i,pi], reason = searchf(dims, optfunc)
      times[i,pi] = toq()
      reason_counts[pi][reason] = get(reason_counts[pi], reason, 0) + 1

      # Print to csv file
      println(csvfh, join([experiment, strftime("%Y%m%d", start_time),
        strftime("%H:%M.%S", start_time), i, prob.name, dims, reason,
        fevals[i,pi], times[i,pi], fbests[i,pi]], ","))
      flush(csvfh)
    end
  end

  for(i in 1:num_problems)
    p = problems[i][2]
    println("\nProblem: ", p.name)
    println("Fitness: ", sumstats(fbests[:,i], (x) -> @sprintf("%.2e", x)))
    println("Time: ", sumstats(times[:,i], format_time))
    println("Num. evals: ", sumstats(fevals[:,i]))
    println("Evals / dim: ", sumstats(fevals[:,i] / BlackBoxOptim.numdims(p), int))
    print("Reasons: "); show(reason_counts[i])
    println("")
  end

  close(csvfh)

  return fbests, fevals, times
end

function compare_params(params, searchf, num_runs = 10)
  num_configs = length(params)
  times = zeros(num_runs, num_configs)
  fbests = zeros(num_runs, num_configs)
  fevals = zeros(num_runs, num_configs)

  for(r in 1:num_runs)
    for(i in 1:num_configs)
      tic()
      xb, fb, fe, reason = searchf(params[i]...)
      times[r,i] = toq()
      fbests[r, i] = fb
      fevals[r, i] = float(fe)
    end
  end

  println("\n\nResults per parameter config")
  println("----------------------------")
  sp = sortperm(mean(fbests, 1)[:])
  for(i in 1:num_configs)
    print(i, ". "); show(params[i]); println("")
    println("Fitness: ", sumstats(fbests[:,i], (x) -> @sprintf("%.6e", x)))
    println("Time: ", sumstats(times[:,i], format_time))
    println("Num. evals: ", sumstats(fevals[:,i]))
    println("")
  end

  return times, fbests, fevals
end

function run_test(f, n; mu = false)
  of = xtransform(n, f)
  tic()
  x, f, fevals, reason = cmsa_es(n, of; covarMatrixSampler = EigenCovarSampler, 
    mu = mu,
    trace = true, known_fmin = 0.0)
  t = toq()
  print("x = ")
  show(x)
  println("\ntermination reason: $(reason)")
  println("fitness = $(f)")
  println("time = $(t)")
  println("fevals = $(fevals)")
  println("fevals/dim = $(fevals/n)")
end
