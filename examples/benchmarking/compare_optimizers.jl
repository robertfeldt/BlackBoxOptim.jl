# Run with: julia --color=yes -L ../../src/BlackBoxOptim.jl compare_optimizers.jl ...

using DataFrames
using BlackBoxOptim
using ArgParse
using CPUTime
using Compat
using Statistics
using Printf: @sprintf

global logfilehandle = nothing
function log(color::Symbol, str)
    printstyled(str, color=color)
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
        raise("Unknown command") # FIXME logfilehandle is not closed
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

ProblemSets = Dict{String,Any}(
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

OptimizerSets = Dict{String,Any}(
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
function increase_runs_per_problem(problemname, numdims)
    k = (problemname, numdims)
    RunsPerProblem[k] = get(RunsPerProblem, k, 0) + 1
end

function fitness_for_opt(family::FunctionBasedProblemFamily, NumDimensions::Int,
                         PopulationSize::Int, MaxFuncEvals::Int,
                         Method::Symbol, TraceMode::Symbol = :verbose)

    problem = instantiate(family, NumDimensions)

    CPUtic()
    res = bboptimize(problem; Method = Method,
        NumDimensions = NumDimensions,
        PopulationSize = PopulationSize,
        TraceMode = TraceMode,
        MaxFuncEvals = MaxFuncEvals
    )

    best_fitness(res), CPUtoc()
end

function latest_git_id(git_repo = joinpath(@__DIR__, "../.."))
    strip(read(`git -C $git_repo log --format="%H" -n 1`, String))
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
                    Method = string(method), PopulationSize = popsize,
                    NumFevals = numfevals,
                    GitID = gitid
            )

            log("\n$(probname), n = $(numdims), optimizer = $(string(method)), run $(rep) of $(NumRepetitions)\n")

            start_time = time()
            ftn, dur = fitness_for_opt(prob, numdims, popsize, ceil(Int, numfevals), method)
            df[:ElapsedTime] = dur
            df[:StartTime] = Libc.strftime("%Y-%m-%d %H:%M:%S", start_time)
            df[:Fitness] = ftn

            push!(dfs, df)
        end
    end

    vcat(dfs)
end

using CSV

function read_benchmark_db(filename)
    if isfile(filename)
        return CSV.read(filename)
    else
        @warn "Benchmark results file $filename not found, returning empty frame"
        return DataFrame()
    end
end

function add_rank_per_group(df, groupcols, rankcol, resultcol)
    by(df, groupcols) do subdf
        orderedsubdf = subdf[sortperm(subdf[:,rankcol]), :]
        orderedsubdf[resultcol] = collect(1:size(orderedsubdf,1))
        return orderedsubdf
    end
end

function list_benchmark_db(db, saveResultCsvFile = nothing)
    numrows = size(db, 1)
    println("Number of runs in db: ", numrows)

    if numrows > 1
        # Find min fitness per problem
        minfitpp = by(db, [:Problem, :NumDims]) do df
            DataFrame(
                MinFitness = minimum(df[:Fitness])
            )
        end

        # Add col with order of magnitude worse than min fitness for each run
        db = join(db, minfitpp, on = [:Problem, :NumDims])
        db[:LogTimesWorseFitness] = log10.(db[:Fitness] ./ db[:MinFitness])

        # Calc median fitness and time per problem and method.
        sumdf = by(db, [:Problem, :NumDims, :Method]) do df
            DataFrame(N = size(df, 1),
                      MedianFitness = median(df[:,:Fitness]),
                      MedianTime = median(df[:,:ElapsedTime]))
        end

        # Rank on median fitness and median time for each problem.
        sumdf = add_rank_per_group(sumdf, [:Problem, :NumDims], :MedianFitness, :RankFitness)
        sumdf = add_rank_per_group(sumdf, [:Problem, :NumDims], :MedianTime, :RankTime)

        # Get number of runs and median magnitude worse per method
        permethod = by(db, [:Method]) do df
            DataFrame(
                NumRuns = size(df, 1),
                MedianLogTimesWorseFitness = round(median(df[:, :LogTimesWorseFitness]), digits=1)
            )
        end

        # and merge with table with mean ranks of fitness and time.
        summarydf = by(sumdf, [:Method]) do df
            DataFrame(
                MeanRank = round(mean(df[:RankFitness]), digits=3),
                Num1sFitness = sum(r -> (r == 1) ? 1 : 0, df[:RankFitness]),
                MeanRankTime = round(mean(df[:RankTime]), digits=3),
                Num1sTime = sum(r -> (r == 1) ? 1 : 0, df[:RankTime]),
            )
        end
        df = join(summarydf, permethod, on = :Method)

        # Now sort and print
        sort!(df, [:MeanRank, :MeanRankTime, :Num1sFitness])
        println(df)
        if saveResultCsvFile !== nothing
            CSV.write(saveResultCsvFile, df)
            println("Results written to file: ", saveResultCsvFile)
        end
    end
end

function save_result_database(db, filename)
    CSV.write(filename, db)
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
    return db
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
            optsel = db[:Method] .== string(optmethod)
            for pd in pset
                probname, numdims, popsize, numfevals = pd
                psel = db[:Problem] .== probname
                dsel = db[:NumDims] .== numdims
                df = db[optsel .& psel .& dsel, :]
                benchfitnesses = convert(Vector{Float64}, df[:Fitness])
                benchtimes = convert(Vector{Float64}, df[:ElapsedTime])
                newfs = Float64[]
                newtimes = Float64[]
                prob = BlackBoxOptim.example_problems[probname]
                for r in 1:nreps
                    runnum += 1
                    log("\n$(probname), n = $(numdims), optimizer = $(string(optmethod)), run $(r) of $(nreps) ($(round(100.0 * runnum / totalruns, digits=2))% of total runs)\n")
                    ftn, dur = fitness_for_opt(prob, numdims, popsize, ceil(Int, numfevals), optmethod)
                    push!(newfs, ftn)
                    push!(newtimes, dur)
                end

                log(:blue, "$(probname)($numdims), $(optmethod):\n")

                ftn_oldmed = median(benchfitnesses)
                ftn_newmed = median(newfs)
                ftn_ord = (ftn_oldmed < ftn_newmed) ? "<" : ">"
                log(:white, "Fitness median $(ftn_oldmed)(old) $(ftn_ord) $(ftn_newmed)(new)\n")
                ftn_pval = pvalue(MannWhitneyUTest(benchfitnesses, newfs))
                if ftn_pval > significancelevel
                    log(:white, "  No statistically significant fitness difference in $nreps repetitions (P-value=$ftn_pval)!\n")
                else
                    log(ftn_ord == "<" ? :red : :green, "  Statistically significant fitness difference in $nreps repetitions (P-value=$ftn_pval)!\n")
                end

                time_oldmed = median(benchtimes)
                time_newmed = median(newtimes)
                time_ord = (time_oldmed < time_newmed) ? "<" : ">"
                log(:white, "Time median $(time_oldmed)(old) $(time_ord) $(time_newmed)(new)\n")
                time_pval = pvalue(MannWhitneyUTest(benchtimes, newtimes))
                if time_pval > significancelevel
                    log(:white, "  No statistically significant time difference in $nreps repetitions (P-value=$time_pval)!\n")
                else
                    log(time_ord == "<" ? :red : :green, "  Statistically significant time difference in $nreps repetitions (P-value=$time_pval)!\n")
                end

                push!(dfs, DataFrame(
                    Problem = probname, NumDims = numdims, Method = string(optmethod), PvalueThreshold = significancelevel,
                    FitnessOldMedian = ftn_oldmed, FitnessOrder = ftn_ord, FitnessNewMedian = ftn_newmed,
                    FitnessPvalue = ftn_pval, FitnessIsSignif = ftn_pval <= significancelevel,
                    TimeOldMedian = time_oldmed, TimeOrder = time_ord, TimeNewMedian = time_newmed,
                    TimePvalue = time_pval, TimeIsSignif = time_pval <= significancelevel)
                )
            end
        end
        df = vcat(dfs...)

        # Use Benjamini-Hochberg to judge which pvalues are significant given we did
        # many comparisons.
        ftn_pvs = convert(Vector{Float64}, df[:FitnessPvalue])
        df[:FitnessSignificantBH001] = benjamini_hochberg(ftn_pvs, 0.01)
        df[:FitnessSignificantBH005] = benjamini_hochberg(ftn_pvs, 0.05)
        df[:FitnessSignificantBH010] = benjamini_hochberg(ftn_pvs, 0.10)

        time_pvs = convert(Vector{Float64}, df[:TimePvalue])
        df[:TimeSignificantBH001] = benjamini_hochberg(time_pvs, 0.01)
        df[:TimeSignificantBH005] = benjamini_hochberg(time_pvs, 0.05)
        df[:TimeSignificantBH010] = benjamini_hochberg(time_pvs, 0.10)

        CSV.write(Libc.strftime("comparison_%Y%m%d_%H%M%S.csv", time()), df)
    else
        df = CSV.read(comparisonfile)
    end
    sort!(df, [:FitnessPvalue])
    report_below_pvalue(df, col_prefix="Fitness", pvalue=1.00)
    report_below_pvalue(df, col_prefix="Fitness", pvalue=0.05)
    report_below_pvalue(df, col_prefix="Fitness", pvalue=0.01)
    report_below_pvalue(df, col_prefix="Time", pvalue=1.00)
    report_below_pvalue(df, col_prefix="Time", pvalue=0.05)
    report_below_pvalue(df, col_prefix="Time", pvalue=0.01)

    # Report (in color) on number of significant differences after Benjamini-Hochberg
    # correction.
    n_ftn_reg = sum(df[:FitnessSignificantBH005] .& (df[:FitnessOrder] .== "<"))
    n_ftn_imp = sum(df[:FitnessSignificantBH005] .& (df[:FitnessOrder] .== ">"))
    printstyled("\n$n_ftn_reg significant fitness regressions at Benjamini-Hochberg 0.05 level\n",
                color=n_ftn_reg > 0 ? :red : :green)
    printstyled("\n$n_ftn_imp significant fitness improvments at Benjamini-Hochberg 0.05 level\n",
                color=n_ftn_imp > 0 ? :green : :white)
    n_time_reg = sum(df[:TimeSignificantBH005] .& (df[:TimeOrder] .== "<"))
    n_time_imp = sum(df[:TimeSignificantBH005] .& (df[:TimeOrder] .== ">"))
    printstyled("\n$n_time_reg significant time regressions at Benjamini-Hochberg 0.05 level\n",
                color=n_time_reg > 0 ? :red : :green)
    printstyled("\n$n_time_imp significant time improvments at Benjamini-Hochberg 0.05 level\n",
                color=n_time_imp > 0 ? :green : :white)
end

function benjamini_hochberg(pvals, alpha = 0.05)
    n = length(pvals)
    (n <= 1) && return pvals # no pvalue correction needed

    perm = sortperm(pvals)
    origperm = invperm(perm)

    thresholds = alpha .* (collect(1:n) ./ n)
    khat = findlast(i -> pvals[perm[i]] < thresholds[i], 1:n)
    (khat === nothing) && return fill(false, n) # There were no significant pvals so return all false

    significant = vcat(fill(true, khat), fill(false, n-khat))
    return significant[origperm]
end

function report_below_pvalue(df; col_prefix="Fitness", pvalue = 0.05)
    selection = df[Symbol(col_prefix, "Pvalue")] .< pvalue
    log("Num problems with $col_prefix p-values < $(pvalue): $(sum(selection))\n")
    # workaround sum(isequal, []) throws
    num_new_worse = sum(df[selection, Symbol(col_prefix, "Order")] .== "<")
    num_new_better = sum(df[selection, Symbol(col_prefix, "Order")] .== ">")
    log_color = num_new_better > num_new_worse ? :green : :red
    log(log_color, "  $col_prefix is better in current run in $(num_new_better) tests\n")
    log(log_color, "  $col_prefix is worse in current run in $(num_new_worse) tests\n")
    any(selection) && println(df[selection, :])
end

@CPUtime main(ARGS)
