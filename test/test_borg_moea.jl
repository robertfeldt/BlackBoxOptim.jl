@testset "BorgMOEA" begin
    @testset "Schaffer1" begin
        res = bboptimize(BlackBoxOptim.Schaffer1Family; Method=:borg_moea,
                         FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                         SearchRange=(-10.0, 10.0), NumDimensions=2, ϵ=0.01,
                         MaxSteps=5000, TraceMode=:silent)
        @test BlackBoxOptim.IGD(BlackBoxOptim.Schaffer1Family.opt_value, pareto_frontier(res),
                                fitness_scheme(res), Val{length(best_candidate(res))}) < 0.05
    end
    @testset "CEC09_UP8" begin
        res = bboptimize(BlackBoxOptim.CEC09_Unconstrained_Set[8]; Method=:borg_moea,
                         NumDimensions=4, ϵ=0.025,
                         MaxSteps=50000, TraceMode=:silent)
        @test BlackBoxOptim.IGD(BlackBoxOptim.CEC09_Unconstrained_Set[8].opt_value, pareto_frontier(res),
                                fitness_scheme(res), Val{length(best_candidate(res))}) < 0.13
    end
    @testset "MaxRestarts option" begin
        res = bboptimize(BlackBoxOptim.Schaffer1Family; Method=:borg_moea,
                         FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                         SearchRange=(-10.0, 10.0), NumDimensions=2, ϵ=0.01,
                         MaxRestarts=3, MaxSteps=100000, TraceMode=:silent)
        @test BlackBoxOptim.general_stop_reason(res) == "Max number of restarts (3) reached"
    end
end
