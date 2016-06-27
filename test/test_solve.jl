facts("Testing solver") do
    function rosenbrock(x)
      return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
    end
    result1 = bboptimize(rosenbrock, SearchRange = (-5.0, 5.0), NumDimensions = 5, Method = :de_rand_1_bin, FitnessTolerance = 1e-5, OptimizationValue = 0.0, MaxSteps = 1000000, TraceMode = :silent);
    result2 = bboptimize(rosenbrock, SearchRange = (-5.0, 5.0), NumDimensions = 5, Method = :de_rand_1_bin, FitnessTolerance = 1e-5, MaxSteps = 1000000, TraceMode = :silent);
    @fact best_fitness(result1) --> less_than(1e-5)
    @fact result1.iterations --> less_than(result2.iterations)
end
