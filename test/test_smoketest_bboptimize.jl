facts("bboptimize() single-objective methods smoketest") do
  rosenbrock2d(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

  for m in keys(BlackBoxOptim.SingleObjectiveMethods)
    context("$(m)") do
      ctrl = bbsetup(rosenbrock2d; Method = m,
        SearchRange = [(-5.0, 5.0), (-2.0, 2.0)], TraceMode = :silent)
      # run first iteration before the main run to exclude compilation from timing
      bboptimize(ctrl, MaxSteps = 1)
      res = bboptimize(ctrl, MaxTime = 0.3)
      @fact length(best_candidate(res)) --> 2
      f = best_fitness(res)
      @fact typeof(f) --> Float64
      @fact f --> less_than(100.0) # this can't be very tight since we give very little time for optimization...
    end
  end
end

facts("bboptimize() multi-objective methods smoketest") do
  schaffer1(x) = (sumabs2(x), sumabs2(x .- 2.0))

  for m in keys(BlackBoxOptim.MultiObjectiveMethods)
    context("$(m)") do
      ctrl = bbsetup(schaffer1; Method = m, MaxSteps=100000,
        SearchRange = [(-10.0, 10.0), (-10.0, 10.0)], TraceMode = :silent,
        FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true), ϵ=0.01)
      # run first iteration before the main run to exclude compilation from timing
      bboptimize(ctrl, MaxSteps = 1)
      res = bboptimize(ctrl, MaxTime = 0.3)
      @fact length(best_candidate(res)) --> 2
      f = best_fitness(res)
      @fact typeof(f) --> NTuple{2,Float64}
      @fact f[1] --> roughly(2.0, atol=5E-2)
      @fact f[2] --> roughly(2.0, atol=5E-2)
    end
  end
end
