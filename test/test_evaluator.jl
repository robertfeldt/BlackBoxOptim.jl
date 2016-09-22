function evaluator_tests(make_eval::Function)
    @testset "Basic evaluation with a single-objective function to be minimized" begin
        e = make_eval()
        @test BlackBoxOptim.num_evals(e) == 0

        fit1 = fitness([0.0, 1.0], e)
        @test BlackBoxOptim.num_evals(e) == 1
        @test BlackBoxOptim.last_fitness(e) == fit1

        fit2 = fitness([2.0, 1.0], e)
        @test BlackBoxOptim.num_evals(e) == 2
        @test BlackBoxOptim.last_fitness(e) == fit2

        @test BlackBoxOptim.numdims(e) == 2

        @test BlackBoxOptim.is_better(0.0, 1.0, e)
        @test BlackBoxOptim.is_better(0.0, -1.0, e) == false

        a, b = [0.0, 0.5], [1.0, 1.0]

        @test BlackBoxOptim.is_better(a, 1.0, e)
        @test BlackBoxOptim.last_fitness(e) == 0.25

        @test BlackBoxOptim.is_better(b, 1.0, e) == false
        @test BlackBoxOptim.last_fitness(e) == 2.0

        @test BlackBoxOptim.best_of(a, b, e) == (a, 0.25)
        @test BlackBoxOptim.best_of(b, a, e) == (a, 0.25)
        BlackBoxOptim.shutdown!(e)
    end

    @testset "update_fitness!()" begin
        e = make_eval()

        candidates = [random_candidate(2, -1.0, 1.0) for i in 1:10]
        BlackBoxOptim.update_fitness!(e, candidates)
        @test BlackBoxOptim.num_evals(e) == length(candidates)
        @test all(c -> isfinite(c.fitness), candidates)
        BlackBoxOptim.shutdown!(e)
    end

    @testset "rank_by_fitness!()" begin
        e = make_eval()

        candidates = [random_candidate(2, -1.0, 1.0) for i in 1:10]
        # partially evaluate fitness
        BlackBoxOptim.update_fitness!(e, candidates[1:5])
        @test BlackBoxOptim.num_evals(e) == 5
        # complete fitness evaluation and sort by it
        BlackBoxOptim.rank_by_fitness!(e, candidates)
        @test BlackBoxOptim.num_evals(e) == 10
        @test sortperm(candidates, by = fitness) == collect(1:10)
        BlackBoxOptim.shutdown!(e)
    end
end

@testset "Evaluator" begin
    # Set up a small example problem
    f(x) = sum(abs2, x)
    p = minimization_problem(f, "", (-1.0, 1.0), 2)
    @testset "ProblemEvaluator" begin
        evaluator_tests(() -> BlackBoxOptim.ProblemEvaluator(p))
    end

    @testset "rank_by_fitness!()" begin
        e = BlackBoxOptim.ProblemEvaluator(p)

        candidates = [random_candidate(2, -1.0, 1.0) for i in 1:10]
        # partially evaluate fitness
        BlackBoxOptim.update_fitness!(e, candidates[1:5])
        @test BlackBoxOptim.num_evals(e) == 5
        # complete fitness evaluation and sort by it
        BlackBoxOptim.rank_by_fitness!(e, candidates)
        @test BlackBoxOptim.num_evals(e) == 10
        @test sortperm(candidates, by = fitness) == collect(1:10)
    end

    @testset "ParallelEvaluator" begin
        using Distributed

        evaluator_tests(() -> BlackBoxOptim.ParallelEvaluator(p, pids=workers()))

        @testset "multi-objective problem" begin
            schaffer1(x) = (sum(abs2, x), sum(xx -> abs2(xx - 2.0), x))
            p = BlackBoxOptim.FunctionBasedProblem(schaffer1, "Schaffer1", ParetoFitnessScheme{2}(is_minimizing=true),
                                                   symmetric_search_space(5, (-10.0, 10.0)))
            a = EpsBoxArchive(EpsBoxDominanceFitnessScheme(fitness_scheme(p)), max_size=100)

            e = BlackBoxOptim.ParallelEvaluator(p, a, pids=workers())
            fit1 = fitness([0.0, 1.0, 2.0, 3.0, 4.0], e)
            @test BlackBoxOptim.num_evals(e) == 1
            @test BlackBoxOptim.last_fitness(e) == fit1.orig

            fit2 = fitness([0.0, -1.0, -2.0, -3.0, -4.0], e)
            @test BlackBoxOptim.num_evals(e) == 2
            @test BlackBoxOptim.last_fitness(e) == fit2.orig

            BlackBoxOptim.shutdown!(e)
        end
    end

end
