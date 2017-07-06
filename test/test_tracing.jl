@testset "Testing methods diagnostic tracing" begin
        rosenbrock(x) = 100.0*sum(i -> abs2(x[i+1] - x[i]^2), 1:length(x)-1) + sum(i -> abs2(x[i] - 1.0), 1:length(x)-1)
        schaffer1(x) = (sum(abs2, x), sum(xx -> abs2(xx - 2.0), x))

        @testset "trace_state()" begin
                for mode in [:silent, :compact, :verbose]
                        @testset "$mode" begin
                                for method in keys(BlackBoxOptim.SingleObjectiveMethods)
                                        opt = bbsetup(rosenbrock; Method=method,
                                                                SearchRange = (-5.0, 5.0), NumDimensions = 2,
                                                                MaxSteps = 1, TraceMode = :silent)
                                        BlackBoxOptim.run!(opt)
                                        BlackBoxOptim.trace_state(DevNull, BlackBoxOptim.optimizer(BlackBoxOptim.lastrun(opt)), mode)
                                end
                                for method in keys(BlackBoxOptim.MultiObjectiveMethods)
                                        opt = bbsetup(schaffer1; Method=method,
                                                                SearchRange = [(-10.0, 10.0), (-10.0, 10.0)],
                                                                FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true), ϵ=0.01,
                                                                MaxSteps = 1, TraceMode = :silent)
                                        BlackBoxOptim.run!(opt)
                                        BlackBoxOptim.trace_state(DevNull, BlackBoxOptim.optimizer(BlackBoxOptim.lastrun(opt)), mode)
                                end
                        end
                end
        end
end
