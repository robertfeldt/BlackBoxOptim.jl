facts("Top-level interface") do
  context("run a simple optimization with mostly defaults") do
    res = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
      MaxTime = 0.5, TraceMode = :silent)
    @fact bestfitness(res) < 0.01 => true
    xbest = best(res)
    @fact typeof(xbest) => Vector{Float64}
    @fact length(xbest) => 2
  end

  context("continue running an optimization after it finished") do
    optctrl = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
      MaxTime = 0.5, TraceMode = :silent)
    res1 = bboptimize(optctrl)
    @fact numruns(optctrl) => 1

    res2 = bboptimize(optctrl; MaxTime = 1.0)
    @fact numruns(optctrl) => 2

    @fact best_fitness(res2) <= best_fitness(res1) => true

    # parameters should be the same except for MaxTime
    for p in keys(parameters(res1))
      if p != :MaxTime
        @fact parameters(res1)[p] == parameters(res2)[p] => true
      end
    end
    @fact parameters(res1)[:MaxTime] => 0.5
    @fact parameters(res2)[:MaxTime] => 1.0
  end

  context("continue running an optimization after serializing to disc") do
    optctrl = bbsetup(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
      MaxTime = 0.5, TraceMode = :silent)
    res1 = bboptimize(optctrl)

    tempfilename = "./temp" * string(rand(1:int(1e8))) * ".tmp"
    open(tempfilename, "w") do fh
      serialize(fh, optctrl)
    end

    # Try to make sure its not in mem:
    opctrl = nothing; gc()

    open(tempfilename, "w") do fh
      ocloaded = deserialize(fh)
    end

    @fact numruns(ocloaded) => 1
    res2 = bboptimize(ocloaded; MaxTime = 1.0)
    @fact numruns(ocloaded) => 2

    @fact bestfitness(res2) <= bestfitness(res1) => true
  end

  context("custom user termination decision") do
    # You can implement your own termination condition by inserting a callback function that returns
    # an OptimizationAction (such as, for example, StopOptimization) or just nothing for no action.
    terminate_func = (optstate) -> numsteps(optstate) >= 100 ? BlackBoxOptim.StopOptimization : nothing
    res = bboptimize(rosenbrock; MaxSteps = 500, SearchRange = (-5.0, 5.0), NumDimensions = 2,
      TraceMode = :silent, StepCallback => terminate_func)
    @fact numsteps(res) => 100
  end

  context("accessing the trace of the optimization") do
    res = bboptimize(rosenbrock; SearchRange = (-5.0, 5.0), NumDimensions = 2,
      MaxTime = 0.5, TraceMode = :silent)
    tr = trace(res)
    numimprovs = length(tr)
    # Last fitness value in trace is the best fitness achieved
    @fact bestfitness(res) == tr[numimprovs] => true
    if numimprovs > 1
      @fact tr[numimprovs] < tr[numimprovs-1] => true
    end
  end
end
