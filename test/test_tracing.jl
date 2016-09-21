@testset "Testing methods diagnostic tracing" begin
        rosenbrock(x) = 100.0*sumabs2(x[2:end] - x[1:end-1].^2) + sumabs2(x[1:end-1]-1.0)
        schaffer1(x) = (sumabs2(x), sumabs2(x .- 2.0))

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
