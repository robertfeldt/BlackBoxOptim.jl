using DataFrames
using BlackBoxOptim
using ArgParse

function main(args)

  s = ArgParseSettings(description = "comparing optimizers to benchmark runs",
    commands_are_required = true,
    version = "0.1.0",
    add_version = true
  )

  @add_arg_table s begin
    "update"
      action = :command
      help = "update benchmark runs"

    "list"
      action = :command
      help = "list info about benchmark runs"

    "compare"
      action = :command
      help = "compare current optimizer(s) to benchmark runs"
  end

  for cmd in ["update", "compare", "list"]
    @add_arg_table s[cmd] begin
      "--benchmarkfile"
        arg_type = String
        default = "benchmark_runs.csv"
        help = "name of benchmark runs db file"
    end
  end

  for cmd in ["update", "compare"]
    @add_arg_table s[cmd] begin  
      "--problems"
        arg_type = String
        default = "all"
        help = "name of problem set"
  
      "--optimizers"
        arg_type = String
        default = "de"
        help = "name of optimizer or optimizer set"
  
      "--numreps"
        arg_type = Integer
        default = 10
        help = "number of repetitions per problem and optimizer"
  
      "--adaptive"
        action = :store_true
        help = "adaptive number of reps until we can prove there is a difference (with a max of 30 reps)"
    end
  end

  pargs = parse_args(args, s)

  cmd = pargs["%COMMAND%"]

  benchmarkfile = pargs[cmd]["benchmarkfile"]

  local pset, optimizers, nreps

  if in(cmd, ["update", "compare"])
    if haskey(ProblemSets, pargs[cmd]["problems"])
      pset = ProblemSets[pargs[cmd]["problems"]]
    else
      ps = pargs[cmd]["problems"]
      throw("Unknown problem set $(ps), existing ones are: ", keys(ProblemSets))
    end

    if haskey(OptimizerSets, pargs[cmd]["optimizers"])
      optimizers = OptimizerSets[pargs[cmd]["optimizers"]]
    else
      os = pargs[cmd]["optimizers"]
      throw("Unknown optimizer set $(os), existing ones are: ", keys(OptimizerSets))
    end

    nreps = pargs[cmd]["numreps"]
  end

  if cmd == "update"
    db = read_benchmark_db(benchmarkfile)
    db = update_benchmarks(db, pset, optimizers, nreps)
    save_result_database(db, benchmarkfile)
  elseif cmd == "list"
    db = read_benchmark_db(benchmarkfile)
    list_benchmark_db(db)
  elseif cmd == "compare"
    db = read_benchmark_db(benchmarkfile)
    compare_optimizers_to_benchmarks(db, pset, optimizers, nreps, 0.01)
  else
    raise("Unknown command")
  end
end

ProblemSets = {
  "easy" => [
    ("Sphere",        5, 20, 5e3),
    ("Sphere",       10, 20, 1e4),
    ("Sphere",       30, 20, 3e4),

    ("Schwefel2.22",  5, 20, 5e3),
    ("Schwefel2.22", 10, 20, 1e4),
    ("Schwefel2.22", 30, 20, 3e4),

    ("Schwefel2.21",  5, 20, 5e3),
    ("Schwefel2.21", 10, 20, 1e4),
    ("Schwefel2.21", 30, 20, 3e4)
  ],

  "harder" => [
    # Harder problems
    ("Schwefel1.2",   5, 20, 5e3),
    ("Schwefel1.2",  10, 50, 5e4),
    ("Schwefel1.2",  30, 50, 2e5),
    ("Schwefel1.2",  50, 50, 3e5),

    ("Rosenbrock",    5, 20, 1e4),
    ("Rosenbrock",   10, 50, 5e4),
    ("Rosenbrock",   30, 50, 2e5),
    ("Rosenbrock",   50, 40, 3e5),
  ],

  "test" => [
    ("Rosenbrock",   30, 50, 2e5),
  ]
}
ProblemSets["all"] = vcat(values(ProblemSets)...)

OptimizerSets = {
  "de" => [:de_rand_1_bin, :de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin, :adaptive_de_rand_1_bin_radiuslimited],
  "test" => [:de_rand_1_bin],
}
OptimizerSets["all"] = collect(keys(BlackBoxOptim.ValidMethods))

function fitness_for_opt(problem, numDimensions, populationSize, numSteps, method)
  problem = BlackBoxOptim.as_fixed_dim_problem(problem, numDimensions)

  println("\n$(problem.name), n = $(numdims(problem)), optimizer = $(string(method))")

  best, fitness = bboptimize(problem; method = method, parameters = {
    :NumDimensions => numDimensions,
    :PopulationSize => populationSize,
    :MaxSteps => numSteps
    })

  fitness
end

function latest_git_id()
  strip(readall(`git log --format="%H" -n 1`))
end

# Test an optimizer on multiple problems and dimensions, return results in dict.
function multitest_opt(problemDescriptions, method; NumRepetitions = 3)
  dfs = DataFrame[]

  gitid = latest_git_id()

  for pd in problemDescriptions
    println("Running $(NumRepetitions) reps for $(method)")
    for rep in 1:NumRepetitions
      probname, numdims, popsize, numsteps = pd
      prob = BlackBoxOptim.example_problems[probname]

      df = DataFrame(Problem = probname, NumDims = numdims,  
        Method = method, PopulationSize = popsize, NumSteps = numsteps, 
        GitID = gitid)

      start_time = time()
      tic()
      ftn = fitness_for_opt(prob, numdims, popsize, numsteps, method)
      df[:ElapsedTime] = toc()
      df[:StartTime] = strftime("%Y-%m-%d %H:%M:%S", start_time)
      df[:Fitness] = ftn

      push!(dfs, df)
    end
  end

  vcat(dfs)
end

function read_benchmark_db(filename)
  if isfile(filename)
    return readtable(filename)
  else
    return DataFrame()
  end
end

function list_benchmark_db(db)
  numrows = size(db, 1)
  println("Number of runs in db: ", numrows)
  if numrows > 1
    sumdf = by(db, [:Problem, :NumDims, :Method], df ->
      DataFrame(N = size(df, 1), MeanFitness = mean(df[:,:Fitness]), 
        StdFitness = std(df[:,:Fitness]), MedianFitness = median(df[:,:Fitness])))

    dfs = DataFrame[]
    for subdf in groupby(sumdf, [:Problem, :NumDims])
      mfitness = subdf[:,:MedianFitness]
      sp = sortperm(mfitness)
      subdf = subdf[sp, :]
      #sort!(subdf, [:MedianFitness])
      n = size(subdf, 1)
      subdf[:Rank] = collect(1:n)
      push!(dfs, subdf)
    end
    sumdf = vcat(dfs...)
    #println(sumdf)

    # Get number of runs per method
    numrunspermethod = by(db, [:Method], df -> DataFrame(NumRuns = size(df, 1)))

    # and merge with table with mean ranks.
    summarydf = by(sumdf, [:Method], df -> DataFrame(MeanRank = mean(df[:,:Rank]), 
      Num1s = sum(map(r -> (r == 1) ? 1 : 0, df[:,:Rank]))))
    df = join(summarydf, numrunspermethod, on = :Method)
    sort!(df; cols = [:MeanRank, :Num1s])
    println(df)
  end
end

function save_result_database(db, filename)
  writetable(filename, db)
end

function update_benchmarks(db, pset, optimizers, nreps = 10)
  for optmethod in optimizers
    @time res = multitest_opt(pset, optmethod; NumRepetitions = nreps)
    if size(db, 1) < 1
      db = res
    else
      db = vcat(db, res)
    end
  end
  db
end

using HypothesisTests

function report_fitness_difference(oldfitnesses, newfitnesses, significancelevel = 0.05)
  mdb, m = median(benchfitnesses), median(fs)
  if mdb > m
    println("  Median fitness of current implementation is LOWER than one in db: ",
      @sprintf("%.3e < %.3e", m, mdb))
  else
    println("  Median fitness of current implementation is HIGHER than one in db: ",
      @sprintf("%.3e < %.3e", m, mdb))
  end
end

function compare_optimizers_to_benchmarks(db, pset, optimizers, nreps, significancelevel::Float64 = 0.05)
  dfs = DataFrame[]
  for optmethod in optimizers
    optsel = db[:,:Method] .== string(optmethod)
    for pd in pset
      probname, numdims, popsize, numsteps = pd
      psel = db[:,:Problem] .== probname
      dsel = db[:,:NumDims] .== numdims
      df = db[optsel & psel & dsel,:]
      benchfitnesses = df[:,:Fitness]
      newfs = Float64[]
      prob = BlackBoxOptim.example_problems[probname]
      for r in 1:nreps
        ftn = fitness_for_opt(prob, numdims, popsize, numsteps, optmethod)
        push!(newfs, ftn)
      end
      pval = pvalue(MannWhitneyUTest(benchfitnesses, newfs))
      println("$(probname)($numdims), $(optmethod):")
      local statsigndiff
      if pval > significancelevel
        println("No statistically significant difference in $nreps repetitions!")
        statsigndiff = "No"
      else
        println("Statistically significant difference in $nreps repetitions!")
        statsigndiff = "Yes"
      end
      oldmed = median(benchfitnesses)
      newmed = median(newfs)
      ord = (oldmed < newmed) ? "<" : ">"
      push!(dfs, DataFrame(Problem = probname, NumDims = numdims, Method = optmethod, 
        OldMedian = oldmed, Order = ord, NewMedian = newmed, 
        Pvalue = pval, Level = significancelevel, StatSignDiff = statsigndiff))
    end
  end
  df = vcat(dfs...)
  println(df)
  writetable(strftime("comparison_%Y%m%d_%H%M%S.csv", time()), df)
end

main(ARGS)