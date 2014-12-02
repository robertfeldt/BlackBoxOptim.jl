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
      method = :de_rand_1_bin, parameters = {:ShowTrace => false, :MaxSteps => 25001})
    @fact fitness < 0.001 => true
  end

  context("Fixed-dimensional problem takes precedence over search range and related params") do
    prob = BlackBoxOptim.fixeddim_problem((x) -> sum(x); range = (10.0, 20.0), dims = 3)
    xbest, rest = bboptimize(prob; search_range = (0.0, 2.0), dimensions = 2, parameters = {:ShowTrace => false})
    @fact length(xbest) => 3
    @fact xbest[1] >= 10.0 => true
    @fact xbest[2] >= 10.0 => true
    @fact xbest[3] >= 10.0 => true
  end

  # Deafult params override this behavior!!! We need to simplify!
  #context("Specify search space only with search_space param") do
  #  xbest, rest = bboptimize((x) -> sum(x); search_space = [(0.0, 2.0), (3.0, 5.0), (10.0, 20.0)], parameters = {:ShowTrace => false})
  #  @fact length(xbest) => 3
  #  @fact (0.0  <= xbest[1] <= 2.0)  => true
  #  @fact (3.0  <= xbest[2] <= 5.0)  => true
  #  @fact (10.0 <= xbest[3] <= 20.0) => true
  #end

  context("fault handling: anydimensional problem") do
    @fact_throws bboptimize(BlackBoxOptim.anydim_problem("dummy", (x) -> sum(x), 0.0:1.0))
  end

  context("fault handling: func & search range but not num dimensions") do
    @fact_throws BlackBoxOptim.setup_problem((x) -> sum(x); search_range = (0.0, 1.0))
  end

  context("restarting an optimizer again and again should gradually improve") do
    optimizer, problem, params = BlackBoxOptim.setup_bboptimize(rosenbrock2d; 
      search_range = (-5.0, 5.0), dimensions = 2, 
      parameters = {:MaxSteps => 10, :ShowTrace => false})
    best10, fitness10, termination_reason10, elapsed_time10, params, num_evals10 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);
    best20, fitness20, termination_reason20, elapsed_time20, params, num_evals20 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);
    params[:MaxSteps] = 980
    best1000, fitness1000, termination_reason1000, elapsed_time1000, params, num_evals1000 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);
    params[:MaxSteps] = 1000
    fitness10000 = best10000 = elapsed_time1000b = 1 # Just so saved outside of loop body...
    for i in 1:9
      best10000, fitness10000, termination_reason10000, elapsed_time1000b, params, num_evals10000 = BlackBoxOptim.run_optimizer_on_problem(optimizer, problem; parameters = params);
    end

    # Fitness is not worse in sub-sequent runs
    @fact (fitness10 >= fitness20 >= fitness1000 >= fitness10000) => true

    # and it should (almost) always be better after 10000 steps than after 10:
    @fact (fitness10 > fitness10000) => true
  end
end