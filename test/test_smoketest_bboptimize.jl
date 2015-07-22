function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

facts("bboptimize smoketest") do
  for(m in keys(BlackBoxOptim.ValidMethods))
    context("testing $(m) method to ensure it works") do
      res = bboptimize(rosenbrock2d; Method = m,
        SearchRange = [(-5.0, 5.0), (-2.0, 2.0)], MaxTime = 0.3, ShowTrace = false)
      @fact length(best_candidate(res)) => 2
      f = best_fitness(res)
      @fact typeof(f) => Float64
      @fact f => less_than(100.0) # this can't be very tight since we give very little time for optimization...
    end
  end
end
