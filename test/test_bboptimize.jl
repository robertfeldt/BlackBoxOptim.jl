function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

function rosenbrock(x)
  return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
end

function sphere(x)
  sum(x.^2)
end

facts("bboptimize") do
  context("example 1 from README") do
    best, fitness = bboptimize(rosenbrock2d, (-5.0, 5.0); dimensions = 2, show_trace = false)
    @fact fitness < 0.001 => true
  end

  context("example 2 from README") do
    best, fitness = bboptimize(rosenbrock2d, [(-5.0, 5.0), (-2.0, 2.0)]; show_trace = false)
    @fact fitness < 0.001 => true
  end

  context("example 3 from README") do
    best, fitness = bboptimize(rosenbrock2d, (-5.0, 5.0); dimensions = 2, method = :de_rand_1_bin, show_trace = false)
    @fact fitness < 0.001 => true
  end

  context("example 4 from README") do
    best, fitness = bboptimize(rosenbrock2d, (-5.0, 5.0); dimensions = 2, 
      method = :random_search, max_time = 10.0, show_trace = false)
    @fact fitness < 0.2 => true
  end

  context("example 5 from README") do
    # dimension 30 is too much for xnes, really, need to change to comparison
    # for same time used or number of fitness evaluations... XNES currently has
    # a huge advantage in using more time and more fitness evals.
    BlackBoxOptim.compare_optimizers(rosenbrock, (-5.0, 5.0); dimensions = 30, max_time = 5.0)
  end

  context("comparing optimizers on more than one problem") do
    BlackBoxOptim.compare_optimizers([(rosenbrock, (-5.0, 5.0)), (sphere, (-5.0, 5.0))]; 
      dimensions = 10, max_time = 5.0)
  end

  context("run one longer example in case there is problem with the reporting in long runs") do
    best, fitness = bboptimize(rosenbrock2d, (-5.0, 5.0); dimensions = 2, 
      method = :de_rand_1_bin, show_trace = false, iterations = 25001)
    @fact fitness < 0.001 => true
  end

 context("test each method option in short runs, just to ensure they work") do
   for(m in keys(BlackBoxOptim.ValidMethods))
     b, f = bboptimize(rosenbrock2d, [(-5.0, 5.0), (-2.0, 2.0)]; show_trace = false, 
       method = m, max_time = 0.4, population_size = 20)
     if (m != :random_search) & (m != :de_rand_2_bin)
       @fact f < 1.0 => true
     end
   end
 end
end