facts("Top-level interface") do
  function rosenbrock(x)
    return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
  end

  context("run a simple optimization with mostly defaults") do
    res = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
      MaxSteps = 1000, ShowTrace = false)
    @fact best_fitness(res) < 0.1 => true
    xbest = best_candidate(res)
    @fact typeof(xbest) => Vector{Float64}
    @fact length(xbest) => 2

    # We also mimic some of the Optim.jl api (although it is a bit strange...)
    @fact f_minimum(res) < 0.1 => true
    @fact minimum(res) == xbest => true
    @fact iteration_converged(res) => true
  end
end
