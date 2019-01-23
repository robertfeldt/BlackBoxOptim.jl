using CPUTime

function compare_optimizers(functionOrProblem, parameters::Parameters = EMPTY_PARAMS;
    Methods = BlackBoxOptim.SingleObjectiveMethodNames, kwargs...)

    parameters = chain(convert(ParamsDict, parameters), kwargs2dict(kwargs))

    results = Any[]
    for m in Methods
        CPUtic()
        res = nothing
        try
            res = bboptimize(functionOrProblem, parameters; Method = m)
        catch e
            @warn "$e when running $m"
            push!(results, (m, e, nothing, missing, missing))
        end
        if res !== nothing
            push!(results, (m, "ok", best_candidate(res), best_fitness(res), CPUtoq()))
        end
    end

    sorted = sort( results, by = (t) -> t[4] )

    if get(parameters, :TraceMode, :compact) != :silent
        println("\n********************************************************************************")
        #println(describe(evaluator))
        for i in 1:length(sorted)
            println("$(i). $(sorted[i][1]): $(sorted[i][2]), fitness = $(sorted[i][4]), time = $(sorted[i][5])")
        end
        println("********************************************************************************\n")
    end

    return sorted
end

function compare_optimizers(problems::Dict{Any, OptimizationProblem},
                            parameters::Parameters = EMPTY_PARAMS;
                            Methods = BlackBoxOptim.SingleObjectiveMethodNames, kwargs...)

    parameters = chain(convert(ParamsDict, parameters), kwargs2dict(kwargs))

    # Lets create an array where we will save how the methods ranks per problem.
    ranks = zeros(length(Methods), length(problems))
    fitnesses = zeros(Float64, length(Methods), length(problems))
    times = zeros(Float64, length(Methods), length(problems))

    problems = collect(problems)

    for i in 1:length(problems)
        name, p = problems[i]
        res = compare_optimizers(p, parameters; Methods = Methods)
        for j in 1:length(res)
            method, best, fitness, elapsedtime = res[j]
            index = findfirst(Methods, method)
            ranks[index, i] = j
            fitnesses[index, i] = fitness
            times[index, i] = elapsedtime
        end
    end

    avg_ranks = round(mean(ranks, dims=2), digits=2)
    avg_fitness = round(mean(fitnesses, dims=2), digits=3)
    avg_times = round(mean(times, dims=2), digits=2)

    perm = sortperm(avg_ranks[:])
    println("\nBy avg rank:")
    for i in 1:length(methods)
        j = perm[i]
        print("\n$(i). $(methods[j]), avg rank = $(avg_ranks[j]), avg fitness = $(avg_fitness[j]), avg time = $(avg_times[j]), ranks = ")
        showcompact(ranks[j,:][:])
    end

    perm = sortperm(avg_fitness[:])
    println("\n\nBy avg fitness:")
    for i in 1:length(methods)
        j = perm[i]
        print("\n$(i). $(methods[j]), avg rank = $(avg_ranks[j]), avg fitness = $(avg_fitness[j]), avg time = $(avg_times[j]), ranks = ")
        showcompact(ranks[j,:][:])
    end

    return ranks, fitnesses
end

"""
Summarize a vector of float values by stating its mean, std dev and median.
"""
function report_on_values(desc, v, lpad = "", rpad = "", digits = 3)
    println("$(lpad)$(desc): $(round(mean(v), sigdigits=digits)) (std. dev = $(round(std(v), sigdigits=digits)), median = $(round(median(v), sigdigits=digits)))")
end

"""
Report the number of times each key was encountered in a count `dict`.

Returns a percentage dict calculated while iterating over the counted items.
"""
function count_dict_report(dict, desc, lpad = "", rpad = "")
    println(desc, ":")
    total = sum(collect(values(dict))) # FIXME collect() should not be required
    pdict = Dict()
    for (r, c) in dict
        pdict[r] = round(100.0*c/total, digits=2)
        println(lpad, r, ": ", c, " (", pdict[r], "%)", rpad)
    end
    pdict
end

"""
Print a report based on a result dict from one set of repeated runs of
an optimization method. Returns the success rate, i.e. number of times the
termination reason was "Within fitness tolerance...".
"""
function report_from_result_dict(statsdict)
    println("Method: $(statsdict[:method])")
    pdict = count_dict_report(statsdict[:reasoncounts], "  Termination reasons", "    ")
    report_on_values("Fitness", statsdict[:fitnesses], "  ")
    report_on_values("Time", statsdict[:times], "  ")
    report_on_values("Num function evals", statsdict[:numevals], "  ")

    success_rate = round(get(pdict, "Within fitness tolerance of optimum", 0.0), digits=3)
    println("  Success rate: ", success_rate, "%\n")
    success_rate
end

function rank_result_dicts_by(result_dicts, byfunc, desc; rev = false,
        descsummary = "mean", digits = 3, rpad = "")

    ranked = BlackBoxOptim.Utils.assign_ranks_within_tolerance(result_dicts; by = byfunc, tolerance = 1e-3, rev = rev)
    println("Ranked by $(descsummary) $(desc):")
    for (rank, rd, value) in ranked
        println("  $(rank). $(rd[:method]), $(round(value, sigdigits=digits))$(rpad)")
    end
    println("")
end

function report_on_methods_results_on_one_problem(problem, result_dicts, numrepeats, max_time, ftol)
    println("********************************************************************************\n")

    println("Problem: $(name(problem)), num dims = $(numdims(problem))")
    println("  Num repeats per method = ", numrepeats)
    println("  Fitness tolerance = ", ftol, " (a run is a success if it reaches to within this value of true optimum)")
    println("  Max time budget per run = ", max_time, " secs\n")

    rank_result_dicts_by(result_dicts, (d) -> d[:success_rate], "success rate (to reach within $(ftol) of optimum)";
            descsummary = "median", rev = true, rpad = "%")
    rank_result_dicts_by(result_dicts, (d) -> median(d[:fitnesses]), "fitness"; descsummary = "median")
    rank_result_dicts_by(result_dicts, (d) -> median(d[:times]), "time (in seconds)";
            descsummary = "median", rpad = " secs")
    rank_result_dicts_by(result_dicts, (d) -> round(Int, median(d[:numevals])), "num function evals"; descsummary = "median")

    for rd in result_dicts
        report_from_result_dict(rd)
    end
end

function repeated_bboptimize(numrepeats, problem, dim, methods, max_time, ftol = 1e-5, extraparams::Parameters = EMPTY_PARAMS)
    fp = instantiate(problem, dim)
    result_dicts = ParamsDict[]

    # Just so they are declared
    ps = best_so_far = nothing

    params = chain(extraparams, ParamsDict(:FitnessTolerance => ftol))

    for m in methods
        ts, fs, nes = zeros(numrepeats), zeros(numrepeats), zeros(Int, numrepeats)
        rcounts = Dict{String,Int}()

        results = map(1:numrepeats) do i
            result = bboptimize(fp, params; MaxTime = max_time, Method = m)

            # Save key data about this run:
            ts[i]  = elapsed_time(result)
            fs[i]  = best_fitness(result)
            nes[i] = f_calls(result)

            # And count only the general stop reasons
            reason = general_stop_reason(result)
            rcounts[reason] = 1 + get(rcounts, reason, 0)
        end

        rdict = ParamsDict(:method => m, :fitnesses => fs, :times => ts, :numevals => nes, :reasoncounts => rcounts)
        rdict[:success_rate] = report_from_result_dict(rdict)
        push!(result_dicts, rdict)
    end

    report_on_methods_results_on_one_problem(fp, result_dicts, numrepeats, max_time, ftol)
end
