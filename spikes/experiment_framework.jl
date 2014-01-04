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

  file_prefix = strftime("$(experiment)_%Y%m%d_%H%M%S", time())

  run_csvfile = join([file_prefix, "_runs.csv"])

  summary_csvfh = 1 # Dummy value just so it exists...
  include_run_csv_header = true # Only the first time...
  add_summary_csv_header = true # Only the first time...

  # We add one run per problem (with runid 0) to ensure everything has been
  # compiled. It is not saved in the data though.
  for(i in 0:num_runs)
    # Random order of running each problem
    problem_indices = shuffle(collect(1:num_problems))
    for(pi in problem_indices)
      num, prob = problems[pi]
      dims = BlackBoxOptim.numdims(prob)
      println("Run $(i) of problem $(name(prob))")
      start_time = time()
      tic()
      params, paramheader, x, fb, nf, reason, archive = searchf(prob)
      t = toq()
      if i == 0
        break
      end
      times[i, pi] = t
      fbests[i, pi] = fb
      fevals[i, pi] = nf
      reason_counts[pi][reason] = get(reason_counts[pi], reason, 0) + 1

      # Save fitness history to a csv file
      save_fitness_history_to_csv_file(archive, run_csvfile; 
        header_prefix = join(["Problem,Dimensions,RunId", paramheader], ","), 
        line_prefix = join(["\"$(name(prob))\",$(dims),$(i)", params], ","),
        include_header = include_run_csv_header)
      include_run_csv_header = false # Only the first round...
      println("Saved fitness history to file: $(run_csvfile)")

      # Create the summary_csv file on the first loop through since we now
      # know the number of parameters
      if add_summary_csv_header
        summary_csvfh = csvfile(["Experiment", "Date", "Time", "RunId", 
          "Problem", "Dimension", paramheader,
          "ElapsedTime", "FuncEvals", "TerminationReason", 
          "Fitness"]; filepath = join([file_prefix, "_summary.csv"]))
        add_summary_csv_header = false
      end

      # Print to summary csv file
      println(summary_csvfh, join([experiment, strftime("%Y-%m-%d", start_time),
        strftime("%T", start_time), i, "\"$(name(prob))\"", dims, params,
        times[i,pi], fevals[i,pi], "\"$(reason)\"", fbests[i,pi]], ","))
      flush(summary_csvfh)
    end
  end

  for(i in 1:num_problems)
    p = problems[i][2]
    println("\nProblem: ", name(p))
    println("Dimensions: ", numdims(p))
    println("Fitness: ", sumstats(fbests[:,i], (x) -> @sprintf("%.2e", x)))
    println("Time: ", sumstats(times[:,i], format_time))
    println("Num. evals: ", sumstats(fevals[:,i]))
    println("Evals / dim: ", sumstats(fevals[:,i] / BlackBoxOptim.numdims(p), round))
    print("Reasons: "); show(reason_counts[i])
    println("")
  end

  close(summary_csvfh)

  println("Finished with experiment $(experiment).")

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
