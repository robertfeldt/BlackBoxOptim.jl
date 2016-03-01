function compare_optimizers(functionOrProblem, parameters::Associative = @compat(Dict{Any,Any}());
  Methods = BlackBoxOptim.MethodNames, kwargs...)

  parameters = convert_and_chain(parameters, kwargs)

  results = Any[]
  for(m in Methods)
    tic()
    res = bboptimize(functionOrProblem, parameters; Method = m)
    push!( results,  (m, best_candidate(res), best_fitness(res), toq()) )
  end

  sorted = sort( results, by = (t) -> t[3] )

  if parameters[:TraceMode] != :silent
    println("\n********************************************************************************")
    #println(describe(evaluator))
    for(i in 1:length(sorted))
      println("$(i). $(sorted[i][1]), fitness = $(sorted[i][3]), time = $(sorted[i][4])")
    end
    println("********************************************************************************\n")
  end

  return sorted

end

function compare_optimizers(problems::Dict{Any, OptimizationProblem},
  parameters::Associative = @compat(Dict{Any,Any}());
  Methods = BlackBoxOptim.MethodNames, kwargs...)

  parameters = convert_and_chain(parameters, kwargs)

  # Lets create an array where we will save how the methods ranks per problem.
  ranks = zeros(length(Methods), length(problems))
  fitnesses = zeros(Float64, length(Methods), length(problems))
  times = zeros(Float64, length(Methods), length(problems))

  problems = collect(problems)

  for i in 1:length(problems)
    name, p = problems[i]
    res = compare_optimizers(p, parameters; Methods = Methods)
    for(j in 1:length(res))
      method, best, fitness, elapsedtime = res[j]
      index = findfirst(Methods, method)
      ranks[index, i] = j
      fitnesses[index, i] = fitness
      times[index, i] = elapsedtime
    end
  end

  avg_ranks = round(mean(ranks, 2), 2)
  avg_fitness = round(mean(fitnesses, 2), 3)
  avg_times = round(mean(times, 2), 2)

  perm = sortperm(avg_ranks[:])
  println("\nBy avg rank:")
  for(i in 1:length(methods))
    j = perm[i]
    print("\n$(i). $(methods[j]), avg rank = $(avg_ranks[j]), avg fitness = $(avg_fitness[j]), avg time = $(avg_times[j]), ranks = ")
    showcompact(ranks[j,:][:])
  end

  perm = sortperm(avg_fitness[:])
  println("\n\nBy avg fitness:")
  for(i in 1:length(methods))
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
  println("$(lpad)$(desc): $(signif(mean(v), digits)) (std. dev = $(signif(std(v), digits)), median = $(signif(median(v), digits)))")
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
    pdict[r] = round(100.0*c/total, 2)
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
  println("  Success rate: ", round(pdict["Within fitness tolerance of optimum"], 3), "%\n")
  pdict["Within fitness tolerance of optimum"]
end

function rank_result_dicts_by(result_dicts, byfunc, desc; rev = false,
  descsummary = "mean", digits = 3, rpad = "")

  ranked = BlackBoxOptim.Utils.assign_ranks_within_tolerance(result_dicts; by = byfunc, tolerance = 1e-3, rev = rev)
  println("Ranked by $(descsummary) $(desc):")
  for (rank, rd, value) in ranked
    println("  $(rank). $(rd[:method]), $(signif(value, digits))$(rpad)")
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
  rank_result_dicts_by(result_dicts, (d) -> int(median(d[:numevals])), "num function evals"; descsummary = "median")

  for rd in result_dicts
    report_from_result_dict(rd)
  end

end

function repeated_bboptimize(numrepeats, problem, dim, methods, max_time, ftol = 1e-5, parameters::Parameters = EMPTY_PARAMS)

  fp = BlackBoxOptim.fixed_dim_problem(problem, dim)
  result_dicts = Dict{Symbol,Any}[]

  # Just so they are declared
  ps = best_so_far = nothing

  params = chain(parameters, @compat Dict{Symbol,Any}(:FitnessTolerance => ftol))

  for m in methods

    ts, fs, nes = zeros(numrepeats), zeros(numrepeats), zeros(Int, numrepeats)
    rcounts = @compat Dict{String,Int}("Within fitness tolerance of optimum" => 0)

    for i in 1:numrepeats
      p = fp # BlackBoxOptim.ShiftedAndBiasedProblem(fp)
      best, fs[i], reason, ts[i], ps, nes[i] = bboptimize(p; max_time = max_time,
        method = m, parameters = params)
      rcounts[reason] = 1 + get(rcounts, reason, 0)
    end

    if best_so_far == nothing
      best_so_far = worst_fitness(ps[:Evaluator])
    end

    # ???
    best_so_far =

    rdict = @compat Dict{Symbol,Any}(:method => m, :fitnesses => fs, :times => ts, :numevals => nes, :reasoncounts => rcounts)
    rdict[:success_rate] = report_from_result_dict(rdict)
    push!(result_dicts, rdict)

  end

  report_on_methods_results_on_one_problem(fp, result_dicts, numrepeats, max_time, ftol)

end
