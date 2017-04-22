# Run with: julia --color=yes -L ../../src/BlackBoxOptim.jl compare_optimizers.jl ...
using DataFrames
using BlackBoxOptim
using ArgParse
using CPUTime
using Compat
using PValueAdjust

global logfilehandle = nothing
function log(color::Symbol, str)
    print_with_color(color, str)
    if logfilehandle != nothing
        print(logfilehandle, str)
        flush(logfilehandle)
    end
end
log(str::AbstractString) = log(:white, str)

function main(args)

    start_time = time()
    global logfilehandle
    logfilehandle = open(Libc.strftime("%Y%m%d_%H%M%S.log", start_time), "a+")

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
      "--benchmarkfile", "-b"
        arg_type = AbstractString
        default = "benchmark_runs.csv"
        help = "name of benchmark runs db file"
    end
  end

  for cmd in ["update", "compare"]
    @add_arg_table s[cmd] begin
      "--problems"
        arg_type = AbstractString
        default = "all"
        help = "name of problem set"

      "--optimizers"
        arg_type = AbstractString
        default = "stable"
        help = "name of optimizer or optimizer set"

      "--numreps", "-n"
        arg_type = Int
        default = 10
        help = "number of repetitions per problem and optimizer"

      "--adaptive"
        action = :store_true
        help = "adaptive number of reps until we can prove there is a difference (with a max of 30 reps)"
    end
  end

  @add_arg_table s["compare"] begin
    "--comparisonfile"
      arg_type = AbstractString
      default = ""
      help = "name of comparison db file"
  end

  @add_arg_table s["list"] begin
    "--outfile", "-o"
      arg_type = AbstractString
      default = ""
      help = "name of csv file where result comparison table is saved"
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
    list_benchmark_db(db, pargs[cmd]["outfile"])
  elseif cmd == "compare"
    compare_optimizers_to_benchmarks(benchmarkfile, pset, optimizers, nreps, 0.01,
        pargs[cmd]["comparisonfile"])
  else
    raise("Unknown command")
  end

  close(logfilehandle)

  println("elapsed time:    ", format_time(time() - start_time))
end

function format_time(t)
  if t < 60.0
    @sprintf("%.1f secs", t)
  elseif t < 3600.0
    m = floor(Int, t/60.0)
    @sprintf("%d min%s ", m, (m>1) ? "s" : "") * format_time(t - 60*m)
  else
    h = floor(Int, t/3600.0)
    @sprintf("%d hour%s ", h, (h>1) ? "s" : "") * format_time(t - 3600*h)
  end
end

ProblemSets = Dict{ASCIIString,Any}(
  "easy" => [
    # ProblemName, NumDims, PopSize, MaxFevals
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

    ("Rastrigin",    50, 50, 5e5),
    ("Rastrigin",   100, 90, 8e5),

    ("Ackley",       50, 50, 5e5),
    ("Ackley",      100, 90, 8e5),

    ("Griewank",     50, 50, 5e5),
    ("Griewank",    100, 90, 8e5),
  ],

  "lowdim" => [
    ("Schwefel1.2",   2, 25, 1e4),
    ("Rosenbrock",    2, 25, 1e4),
    ("Rastrigin",     2, 25, 1e4),
    ("Ackley",        2, 25, 1e4),
    ("Griewank",      2, 25, 1e4),
  ],

  "test" => [
    ("Rosenbrock",   30, 50, 2e5),
  ]
)
ProblemSets["all"] = vcat(ProblemSets["easy"], ProblemSets["harder"])

OptimizerSets = Dict{ASCIIString,Any}(
  "de" => [:de_rand_1_bin, :de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin, :adaptive_de_rand_1_bin_radiuslimited],
  "stable_non_de" => [:probabilistic_descent, :generating_set_search, :random_search],
  "nes" => [:xnes, :separable_nes, :dxnes],
  "test" => [:de_rand_1_bin],
)
OptimizerSets["all"] = unique(collect(keys(BlackBoxOptim.SingleObjectiveMethods)))
OptimizerSets["stable"] = vcat(OptimizerSets["de"], OptimizerSets["stable_non_de"])

# Keep track of how many times we have executed each problem for this run of the script,
# since we should never take time measurements before we have ensured everything has been compiled.
const RunsPerProblem = Dict{Any,Int}()
runs_per_problem(problemname, numdims) = get(RunsPerProblem, (problemname, numdims), 0)
increase_runs_per_problem(problemname, numdims) = begin
  k = (problemname, numdims)
  RunsPerProblem[k] = get(RunsPerProblem, k, 0) + 1
end

function fitness_for_opt(family::FunctionBasedProblemFamily, NumDimensions::Int, PopulationSize::Int, MaxFuncEvals::Int,
    Method::Symbol, TraceMode::Symbol = :verbose)

  problem = instantiate(family, NumDimensions)

  res = bboptimize(problem; Method = Method,
    NumDimensions = NumDimensions,
    PopulationSize = PopulationSize,
    TraceMode = TraceMode,
    MaxFuncEvals = MaxFuncEvals
  )

  best_fitness(res)
end

function latest_git_id()
  strip(readall(`git -C "../.." log --format="%H" -n 1`))
end

# Test an optimizer on multiple problems and dimensions, return results in dict.
function multitest_opt(problemDescriptions, method; NumRepetitions = 3)
  dfs = DataFrame[]

  gitid = latest_git_id()

  for pd in problemDescriptions
    log("Running $(NumRepetitions) reps for $(method)\n")
    for rep in 1:NumRepetitions
      probname, numdims, popsize, numfevals = pd

      prob = BlackBoxOptim.example_problems[probname]

      # Ensure everything is compiled before we start measuring the 1st time.
      if runs_per_problem(probname, numdims) == 0
        fitness_for_opt(prob, numdims, popsize, 100, method, :silent)
        increase_runs_per_problem(probname, numdims)
      end

      df = DataFrame(Problem = probname, NumDims = numdims,
        Method = string(method), PopulationSize = popsize, NumFevals = numfevals,
        GitID = gitid)

      log("\n$(probname), n = $(numdims), optimizer = $(string(method)), run $(rep) of $(NumRepetitions)\n")

      start_time = time()
      CPUtic()
      ftn = fitness_for_opt(prob, numdims, popsize, ceil(Int, numfevals), method)
      df[:ElapsedTime] = CPUtoc()
      df[:StartTime] = Libc.strftime("%Y-%m-%d %H:%M:%S", start_time)
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

function add_rank_per_group(df, groupcols, rankcol, resultcol)
    dfs = DataFrame[]
    for subdf in groupby(df, groupcols)
        orderedsubdf = subdf[sortperm(subdf[:,rankcol]), :]
        orderedsubdf[resultcol] = collect(1:size(orderedsubdf,1))
        push!(dfs, orderedsubdf)
    end
    vcat(dfs...)
end

function list_benchmark_db(db, saveResultCsvFile = nothing)
  numrows = size(db, 1)
  println("Number of runs in db: ", numrows)

  if numrows > 1
    # Find min fitness per problem
    minfitpp = by(db, [:Problem, :NumDims], df -> DataFrame(
        MinFitness = minimum(df[:,:Fitness])))

    # Add col with order of magnitude worse than min fitness for each run
    db = join(db, minfitpp, on = [:Problem, :NumDims])
    db[:LogTimesWorseFitness] = log10(db[:Fitness] ./ db[:MinFitness])

    # Calc median fitness and time per problem and method.
    sumdf = by(db, [:Problem, :NumDims, :Method], df ->
        DataFrame(N = size(df, 1),
            MedianFitness = median(df[:,:Fitness]),
            MedianTime = median(df[:,:ElapsedTime])))

    # Rank on median fitness and median time for each problem.
    sumdf = add_rank_per_group(sumdf, [:Problem, :NumDims], :MedianFitness, :RankFitness)
    sumdf = add_rank_per_group(sumdf, [:Problem, :NumDims], :MedianTime, :RankTime)

    # Get number of runs and median magnitude worse per method
    permethod = by(db, [:Method], df -> DataFrame(
        NumRuns = size(df, 1),
        MedianLogTimesWorseFitness = round(median(df[:, :LogTimesWorseFitness]), 1),
        )
    )

    # and merge with table with mean ranks of fitness and time.
    summarydf = by(sumdf, [:Method], df -> DataFrame(
        MeanRank = round(mean(df[:,:RankFitness]), 3),
        Num1sFitness = sum(map(r -> (r == 1) ? 1 : 0, df[:,:RankFitness])),
        MeanRankTime = round(mean(df[:,:RankTime]), 3),
        Num1sTime = sum(map(r -> (r == 1) ? 1 : 0, df[:,:RankTime])),
        )
    )
    df = join(summarydf, permethod, on = :Method)

    # Now sort and print
    sort!(df; cols = [:MeanRank, :MeanRankTime, :Num1sFitness])
    println(df)
    if typeof(saveResultCsvFile) <: AbstractString && length(saveResultCsvFile) > 0
        writetable(saveResultCsvFile, df)
        println("Results written to file: ", saveResultCsvFile)
    end
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
    #save_result_database(db, "temp.csv")
    list_benchmark_db(db)
  end
  db
end

using HypothesisTests

function compare_optimizers_to_benchmarks(benchmarkfile, pset, optimizers, nreps,
    significancelevel::Float64 = 0.05, comparisonfile = "")
    if comparisonfile == ""
        db = read_benchmark_db(benchmarkfile)
        dfs = DataFrame[]
        totalruns = length(optimizers) * nreps * length(pset)
        runnum = 0
        for optmethod in optimizers
            optsel = db[:,:Method] .== string(optmethod)
            for pd in pset
                probname, numdims, popsize, numfevals = pd
                psel = db[:,:Problem] .== probname
                dsel = db[:,:NumDims] .== numdims
                df = db[optsel & psel & dsel,:]
                benchfitnesses = df[:,:Fitness]
                newfs = Float64[]
                prob = BlackBoxOptim.example_problems[probname]
                for r in 1:nreps
                    runnum += 1
                    log("\n$(probname), n = $(numdims), optimizer = $(string(optmethod)), run $(r) of $(nreps) ($(round(100.0 * runnum / totalruns, 2))% of total runs)\n")
                    ftn = fitness_for_opt(prob, numdims, popsize, ceil(Int, numfevals), optmethod)
                    push!(newfs, ftn)
                end
                pval = pvalue(MannWhitneyUTest(benchfitnesses, newfs))
                log(:blue, "$(probname)($numdims), $(optmethod):\n")
                local statsigndiff
                if pval > significancelevel
                    log(:green, "  No statistically significant difference in $nreps repetitions!\n")
                    statsigndiff = "No"
                else
                    log(:red, "  Statistically significant difference in $nreps repetitions!\n")
                    statsigndiff = "Yes"
                end
                oldmed = median(benchfitnesses)
                newmed = median(newfs)
                ord = (oldmed < newmed) ? "<" : ">"
                push!(dfs, DataFrame(Problem = probname, NumDims = numdims, Method = string(optmethod),
                OldMedian = oldmed, Order = ord, NewMedian = newmed,
                Pvalue = pval, Level = significancelevel, StatSignDiff = statsigndiff))
                log(:white, "  (Old median) $(oldmed) $(ord) $(newmed) (New median)\n")
            end
        end
        df = vcat(dfs...)

        # Use Benjamini-Hochberg to judge which pvalues are significant given we did
        # many comparisons.
        pvs = convert(Array, df[:Pvalue])
        @show pvs
        df[:SignificantBH001] = benjamini_hochberg(pvs, 0.01)
        df[:SignificantBH005] = benjamini_hochberg(pvs, 0.05)
        df[:SignificantBH010] = benjamini_hochberg(pvs, 0.10)

        writetable(Libc.strftime("comparison_%Y%m%d_%H%M%S.csv", time()), df)
    else
        df = readtable(comparisonfile)
    end
    sort!(df; cols = [:Pvalue])
    report_below_pvalue(df, 1.00)
    report_below_pvalue(df, 0.05)
    report_below_pvalue(df, 0.01)

    # Report (in color) on number of significant differences after Benjamini-Hochberg
    # correction.
    color = (sum(df[:SignificantBH005]) > 0) ? :red : :green
    print_with_color(color, 
      "\nNum significant at Benjamini-Hochberg 0.05 level: ", string(sum(df[:SignificantBH005])),
      "\n")
end

function benjamini_hochberg(pvals, alpha = 0.05)
  n = length(pvals)
  if n <= 1
    return pvals
  end
  perm = sortperm(pvals)
  origperm = sortperm(perm)
  sortedps = pvals[perm]
  
  tresholds = alpha * collect(1:n) / n
  khat = n+1 - findfirst(reverse(sortedps .<= tresholds))
  if khat > n
    return (ones(n).>10.0) # There were no significant ones so return all false
  else
    significant = vcat(ones(khat), zeros(n-khat)) .> 0.0
    significant[origperm]
  end
end

function report_below_pvalue(df, pvalue = 0.05)
    selection = df[:,:Pvalue] .< pvalue
    if pvalue >= 1.0
        log("Num problems (any p-value) = $(sum(selection))\n")
    else
        log("Num problems (with p-values < $(pvalue)) = $(sum(selection))\n")
    end
    num_new_worse = sum(df[selection, :Order] .== "<")
    num_new_better = sum(df[selection, :Order] .== ">")
    log(:green, "  Num where new implementation is better = $(num_new_better)\n")
    log(:red, "  Num where new implementation is worse = $(num_new_worse)\n")
    if sum(selection) > 0 && pvalue < 1.0
        println(df[selection,:])
    end
end

@CPUtime main(ARGS)
