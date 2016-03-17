facts("BorgMOEA") do
    context("Schaffer1") do
      res = bboptimize(BlackBoxOptim.Schaffer1Family; Method=:borg_moea,
                       FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                       SearchRange=(-10.0, 10.0), NumDimensions=2, ϵ=0.01,
                       MaxSteps=5000, TraceMode=:silent)
      @fact BlackBoxOptim.IGD(BlackBoxOptim.Schaffer1Family.opt_value, pareto_frontier(res),
                fitness_scheme(res), Val{length(best_candidate(res))}) --> less_than(0.05)
    end
    context("CEC09_UP8") do
      res = bboptimize(BlackBoxOptim.CEC09_Unconstrained_Set[8]; Method=:borg_moea,
                       NumDimensions=4, ϵ=0.01,
                       MaxSteps=15000, TraceMode=:silent)
      @fact BlackBoxOptim.IGD(BlackBoxOptim.CEC09_Unconstrained_Set[8].opt_value, pareto_frontier(res),
                fitness_scheme(res), Val{length(best_candidate(res))}) --> less_than(0.13)
    end
end
