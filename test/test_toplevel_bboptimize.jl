facts("Top-level interface") do
  function rosenbrock(x)
    return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
  end

  context("run a simple optimization") do
    context("using bboptimize() with mostly defaults") do
      res = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
        MaxSteps = 2000, TraceMode = :silent)
      @fact best_fitness(res) --> less_than(0.1)
      xbest = best_candidate(res)
      @fact typeof(xbest) --> Vector{Float64}
      @fact length(xbest) --> 2

      # We also mimic some of the Optim.jl api (although it is a bit strange...)
      @fact f_minimum(res) --> less_than(0.1)
      @fact minimum(res) --> xbest
      @fact iteration_converged(res) --> true
    end

    context("using bbsetup()/bboptimize() with mostly defaults") do
      opt = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
        MaxSteps = 2000, TraceMode = :silent)
      @fact numruns(opt) --> 0
      @fact isa(problem(opt), BlackBoxOptim.FunctionBasedProblem) --> true
      res = bboptimize(opt)
      @fact numruns(opt) --> 1
      @fact isa(lastrun(opt), BlackBoxOptim.OptRunController) --> true

      @fact best_fitness(res) --> less_than(0.1)
      xbest = best_candidate(res)
    end

    context("using non-population optimizer") do
      res = bboptimize(rosenbrock; Method=:generating_set_search,
                       SearchRange = (-5.0, 5.0), NumDimensions = 2,
                       MaxSteps = 5000, TraceMode = :silent)
      @fact best_fitness(res) --> less_than(1.0)
      xbest = best_candidate(res)
      @fact typeof(xbest) --> Vector{Float64}
    end

    context("using population optimizer") do
      res = bboptimize(rosenbrock; Method=:adaptive_de_rand_1_bin,
                       SearchRange = (-5.0, 5.0), NumDimensions = 2,
                       MaxSteps = 5000, TraceMode = :silent)
      @fact best_fitness(res) --> less_than(0.1)
      xbest = best_candidate(res)
      @fact typeof(xbest) --> Vector{Float64}
      xpop = population(res)
      @fact isa(xpop, BlackBoxOptim.Population) --> true
      @fact popsize(xpop) --> greater_than(0)
      @fact numdims(xpop) --> 2

      # Access a few individual solution vectors in the population...
      @fact isa(xpop[1], Array{Float64, 1}) --> true # first solution vector
      @fact isa(xpop[popsize(xpop)], Array{Float64, 1}) --> true # first solution vector
      rand_solution_idx = rand(1:popsize(xpop))
      @fact isa(xpop[rand_solution_idx], Array{Float64, 1}) --> true # random solution vector

      # and to access their fitness values:
      @fact isa(fitness(xpop, 1), Float64) --> true
      @fact isa(fitness(xpop, popsize(xpop)), Float64) --> true
      @fact isa(fitness(xpop, rand_solution_idx), Float64) --> true

      # Ensure the lowest fitness value is the one returned by best_fitness
      min_fitness_value = minimum(map(i -> fitness(xpop, i), 1:popsize(xpop)))
      @fact (min_fitness_value == best_fitness(res)) --> true
    end

if BlackBoxOptim.enable_parallel_methods
    context("using population optimizer and parallel evaluator") do
      opt = bbsetup(rosenbrock; Method=:adaptive_de_rand_1_bin,
                    SearchRange = (-5.0, 5.0), NumDimensions = 2,
                    MaxSteps = 2000, TraceMode = :silent, Workers=workers())
      res = bboptimize(opt)
      @fact isa(BlackBoxOptim.evaluator(lastrun(opt)), BlackBoxOptim.ParallelEvaluator) --> true
    end
end

  end

  context("continue running an optimization after it finished") do
    optctrl = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
      MaxTime = 0.5, TraceMode = :silent)

    res1 = bboptimize(optctrl)
    @fact numruns(optctrl) --> 1

    res2 = bboptimize(optctrl; MaxTime = 1.0)
    @fact numruns(optctrl) --> 2

    @fact best_fitness(res2) --> less_than_or_equal(best_fitness(res1))

    # parameters should be the same except for MaxTime
    for p in keys(flatten(parameters(res1)))
      if p != :MaxTime
        @fact parameters(res1)[p] --> parameters(res2)[p]
      end
    end
    @fact parameters(res1)[:MaxTime] --> 0.5
    @fact parameters(res2)[:MaxTime] --> 1.0
  end

  context("return results after interruption") do
      i = 0
      function rosenbrock_throwing(x)
          i += 1
          if i < 50
              return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
          else
              throw(InterruptException())
          end
      end
      context(":RecoverResults on") do
        i = 0
        optctrl = bbsetup(rosenbrock_throwing; SearchRange = (-5.0, 5.0), NumDimensions = 100,
            MaxSteps=100, TraceMode=:silent, RecoverResults=true)
        res = bboptimize(optctrl)
        @fact BlackBoxOptim.stop_reason(res) --> (@sprintf "%s" InterruptException())
      end

      context(":RecoverResults off") do
        i = 0 # reset the counter, otherwise it will throw in the setup
        optctrl = bbsetup(rosenbrock_throwing; SearchRange = (-5.0, 5.0), NumDimensions = 100,
            MaxSteps=100, TraceMode=:silent, RecoverResults=false)
        @fact_throws InterruptException bboptimize(optctrl)
      end
  end

  context("continue running an optimization after serializing to disc") do
    optctrl = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 100,
      MaxTime = 0.5, TraceMode = :silent)
    res1 = bboptimize(optctrl)

    local tempfilename = "./temp" * string(rand(1:10^8)) * ".tmp"

    try # To ensure we delete the temp file afterwards...
      open(tempfilename, "w") do fh
        serialize(fh, optctrl)
      end

      # Try to make sure its not in mem:
      opctrl = nothing; gc()

      local ocloaded
      open(tempfilename, "r") do fh
        ocloaded = deserialize(fh)
      end

      @fact numruns(ocloaded) --> 1
      res2 = bboptimize(ocloaded; MaxTime = 1.0)
      @fact numruns(ocloaded) --> 2

      @fact best_fitness(res2) --> less_than_or_equal(best_fitness(res1))
    finally
      if isfile(tempfilename)
        rm(tempfilename)
      end
    end
  end

  context("TargetFitness option works") do
    # FIXME use the same (fixed?) random seed to guarantee reproducibility
    result1 = bboptimize(rosenbrock, SearchRange = (-5.0, 5.0), NumDimensions = 5,
                         Method = :de_rand_1_bin, FitnessTolerance = 1e-5,
                         MaxSteps = 1000000, TraceMode = :silent,
                         TargetFitness = 0.0)
    result2 = bboptimize(rosenbrock, SearchRange = (-5.0, 5.0), NumDimensions = 5,
                         Method = :de_rand_1_bin, FitnessTolerance = 1e-5,
                         MaxSteps = 1000000, TraceMode = :silent)
    @fact best_fitness(result1) --> less_than(1e-5)
    @fact result1.iterations --> less_than(result2.iterations)
  end
end
