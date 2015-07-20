facts("Top-level interface") do
  context("run a simple optimization with mostly defaults") do
    res = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
      MaxTime = 0.5, ShowTrace = false)
    @fact bestfitness(res) < 0.01 => true
    xbest = best(res)
    @fact typeof(xbest) => Vector{Float64}
    @fact length(xbest) => 2
  end

  context("continue running an optimization after it finished") do
    res1 = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
      MaxTime = 0.5, ShowTrace = false)
    res2 = bboptimize(res1; MaxTime = 1.0) # or should it be bboptimize(optcontroller(res1)...) ?

    @fact bestfitness(res2) <= bestfitness(res1) => true
    @fact numruns(res1) => 1
    @fact numruns(res2) => 2

    # parameters should be the same except for MaxTime
    for p in keys(parameters(res1))
      if p != :MaxTime
        @fact parameters(res1)[p] == parameters(res2)[p] => true
      end
    end
    @fact parameters(res1)[:MaxTime] => 0.5
    @fact parameters(res1)[:MaxTime] => 1.0
  end

  context("custom user termination decision") do
    terminate_func = (optctrl) -> numsteps(optctrl) >= 100 ? stop_optimization!(optctrl) : nothing
    res = bboptimize(rosenbrock; MaxSteps = 500, SearchRange = (-5.0, 5.0), NumDimensions = 2,
      ShowTrace = false, StepCallback => terminate_func)
    @fact numsteps(res) => 100
  end
end
