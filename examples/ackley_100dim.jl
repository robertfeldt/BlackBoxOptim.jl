using BlackBoxOptim

# Summarize a vector of float values by stating its mean, std dev and median.
function report(desc, v, lpad = "", rpad = "", digits = 3)
  println("$(lpad)$(desc): $(signif(mean(v), digits)) (std. dev = $(signif(std(v), digits)), median = $(signif(median(v), digits)))")
end

# Report on the number of times each key in a count dict was encountered. 
# Returns a percentage dict calculated while iterating over the counted items.
function count_dict_report(dict, desc, lpad = "", rpad = "")
  println(desc, ":")
  total = sum(collect(values(dict)))
  pdict = Dict()
  for (r, c) in dict
    pdict[r] = round(100.0*c/total, 2)
    println(lpad, r, ": ", c, " (", pdict[r], "%)", rpad)
  end
  pdict
end

# Print a report based on a result dict from one set of repeated runs of
# an optimization method. Returns the success rate, i.e. number of times the
# termination reason was "Within fitness tolerance...".
function report_from_result_dict(statsdict)
  println("Method: $(statsdict[:method])")
  pdict = count_dict_report(statsdict[:reasoncounts], "  Termination reasons", "    ")
  report("Fitness", statsdict[:fitnesses], "  ")
  report("Time", statsdict[:times], "  ")
  report("Num function evals", statsdict[:numevals], "  ")
  println("  Success rate: ", round(pdict["Within fitness tolerance of optimum"], 3), "%\n")
  pdict["Within fitness tolerance of optimum"]
end

# Assign ranks to values but keep the rank the same if the values are within
# tolerance of each other.
function assign_rank_within(values; byfunc = (x) -> x, tolerance = 1e-5, rev = false)

  perm = sortperm(values, by = byfunc, rev = rev)
  ranked = Any[]
  rank = 1
  prev = byfunc(values[perm[1]])
  num_of_this_rank = 0

  for i in 1:length(perm)
    r = values[perm[i]]
    v = byfunc(r)
    if abs(prev - v) > tolerance
      rank += num_of_this_rank
      num_of_this_rank = 1
    else
      num_of_this_rank += 1
    end
    push!(ranked, (rank, r, v))
    prev = v
  end

  ranked

end

function rank_result_dicts_by(result_dicts, byfunc, desc; rev = false, 
  descsummary = "mean", digits = 3, rpad = "")

  ranked = assign_rank_within(result_dicts; byfunc = byfunc, tolerance = 1e-3, rev = rev)
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

function repeated_bboptimize(numrepeats, problem, dim, methods, max_time, ftol = 1e-5)

  fp = BlackBoxOptim.as_fixed_dim_problem(problem, dim)
  result_dicts = Any[]

  for m in methods

    ts, fs, nes = zeros(numrepeats), zeros(numrepeats), zeros(Int64, numrepeats)
    rcounts = {"Within fitness tolerance of optimum" => 0}

    for i in 1:numrepeats
      p = fp # BlackBoxOptim.ShiftedAndBiasedProblem(fp)
      best, fs[i], reason, ts[i], ps, nes[i] = bboptimize(p; max_time = max_time, 
        method = m, parameters = {:FitnessTolerance => ftol})
      rcounts[reason] = 1 + get(rcounts, reason, 0)
    end

    rdict = {:method => m, :fitnesses => fs, :times => ts, :numevals => nes, :reasoncounts => rcounts}
    rdict[:success_rate] = report_from_result_dict(rdict)
    push!(result_dicts, rdict)

  end

  report_on_methods_results_on_one_problem(fp, result_dicts, numrepeats, max_time, ftol)

end

p = BlackBoxOptim.example_problems["Ackley"]
repeated_bboptimize(5, p, 100, [
  :generating_set_search, 
  :adaptive_de_rand_1_bin_radiuslimited,
  :random_search,
  :separable_nes], 
  20.0)