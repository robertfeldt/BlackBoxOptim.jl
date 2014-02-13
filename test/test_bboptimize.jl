function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

function rosenbrock(x)
  return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
end

facts("bboptimize") do
  context("example 1 from README") do
    best, fitness = bboptimize(rosenbrock2d; search_range = (-5.0, 5.0), dimensions = 2, parameters = {:ShowTrace => false})
    @fact fitness < 0.001 => true
  end

  context("example 2 from README") do
    best, fitness = bboptimize(rosenbrock2d; search_range = [(-5.0, 5.0), (-2.0, 2.0)], parameters = {:ShowTrace => false})
    @fact fitness < 0.001 => true
  end

  context("example 3 from README") do
    best, fitness = bboptimize(rosenbrock2d; search_range = (-5.0, 5.0), dimensions = 2, method = :de_rand_1_bin, parameters = {:ShowTrace => false})
    @fact fitness < 0.001 => true
  end

  context("example 4 from README") do
    best, fitness = bboptimize(rosenbrock2d; search_range = (-5.0, 5.0), dimensions = 2, 
      method = :random_search, max_time = 10.0, parameters = {:ShowTrace => false})
    @fact fitness < 0.2 => true
  end

  context("example 5 from README") do
    BlackBoxOptim.compare_optimizers(rosenbrock; search_range = (-5.0, 5.0), dimensions = 30, max_time = 5.0, parameters = {:ShowTrace => false})
  end

  context("run one longer example in case there is problem with the reporting in long runs") do
    best, fitness = bboptimize(rosenbrock2d; search_range = (-5.0, 5.0), dimensions = 2, 
      method = :de_rand_1_bin, parameters = {:ShowTrace => true, :MaxSteps => 25001})
    @fact fitness < 0.001 => true
  end
end