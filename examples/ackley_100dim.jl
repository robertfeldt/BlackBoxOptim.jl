using BlackBoxOptim
report(desc, v) = println("$(desc) = $(round(mean(v), 4)) +/- $(round(std(v), 4))")
function repeated_bboptimize(numrepeats, problem, dim, method, max_time)
  ts, fs, nes = zeros(numrepeats), zeros(numrepeats), zeros(Int64, numrepeats)
  fp = BlackBoxOptim.as_fixed_dim_problem(problem, dim)
  for i in 1:numrepeats
    p = BlackBoxOptim.ShiftedAndBiasedProblem(fp)
    best, fs[i], termination_reason, ts[i], ps, nes[i] = bboptimize(p; max_time = max_time, method = method)
  end
  println("Method: $(method)")
  report("Fitness", fs)
  report("Time", ts)
  report("Num function evals", nes)
end
p = BlackBoxOptim.example_problems["Ackley"]
#repeated_bboptimize(30, p, 100, :generating_set_search, 30.0)
repeated_bboptimize(30, p, 100, :adaptive_de_rand_1_bin_radiuslimited, 30.0)
