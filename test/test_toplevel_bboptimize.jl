facts("Top-level interface") do
  function rosenbrock(x)
    return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
  end

  context("run a simple optimization with mostly defaults") do
    res = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
      MaxSteps = 2000, ShowTrace = false)
    @fact best_fitness(res) < 0.1 => true
    xbest = best_candidate(res)
    @fact typeof(xbest) => Vector{Float64}
    @fact length(xbest) => 2

    # We also mimic some of the Optim.jl api (although it is a bit strange...)
    @fact f_minimum(res) < 0.1 => true
    @fact minimum(res) == xbest => true
    @fact iteration_converged(res) => true
  end

  context("continue running an optimization after it finished") do
    optctrl = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
      MaxTime = 0.5, ShowTrace = true)

    res1 = bboptimize(optctrl)
    @fact numruns(optctrl) => 1

    res2 = bboptimize(optctrl; MaxTime = 1.0)
    @fact numruns(optctrl) => 2

    @fact best_fitness(res2) <= best_fitness(res1) => true

    # parameters should be the same except for MaxTime
    for p in keys(flatten(parameters(res1)))
      if p != :MaxTime
        @fact parameters(res1)[p] == parameters(res2)[p] => true
      end
    end
    @fact parameters(res1)[:MaxTime] => 0.5
    @fact parameters(res2)[:MaxTime] => 1.0
  end
end
