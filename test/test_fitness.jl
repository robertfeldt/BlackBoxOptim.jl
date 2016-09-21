@testset "Fitness" begin
    @testset "hat_compare() Float64" begin
        @test hat_compare(1.0, 2.0) == -1
        @test hat_compare(-1.0, 1.0) == -1

        @test hat_compare(2.0, 1.0) == 1
        @test hat_compare(-2.0, -3.0) == 1

        @test hat_compare(1.0, 1.0) == 0
        @test hat_compare(0.0, 0.0) == 0
        @test hat_compare(-1.0, -1.0) == 0
    end

    @testset "ScalarFitnessScheme" begin
            @testset "is_minimizing()" begin
                mins = MinimizingFitnessScheme
                @test is_minimizing(mins)

                maxs = MaximizingFitnessScheme
                @test is_minimizing(maxs) == false
            end

            @testset "hat_compare(..., MinimizingFitnessScheme)" begin
                scheme = MinimizingFitnessScheme

                @test hat_compare(1.0, 2.0, scheme) == -1
                @test hat_compare(2.0, 1.0, scheme) == 1
                @test hat_compare(1.0, 1.0, scheme) == 0
            end

            @testset "hat_compare(..., MaximizingFitnessScheme)" begin
                scheme = MaximizingFitnessScheme

                @test hat_compare(1.0, 2.0, scheme) == 1
                @test hat_compare(2.0, 1.0, scheme) == -1
                @test hat_compare(1.0, 1.0, scheme) == 0
            end

            # FIXME enable tests once v0.5 issue #14919 is fixed
            @testset "fitness_scheme(x, y)" begin
                mins = MinimizingFitnessScheme
                #@test !mins(5.0, 3.0)
                #@test mins(1.0, 2.0)

                maxs = MaximizingFitnessScheme
                #@test !maxs(3.0, 5.0)
                #@test maxs(2.0, 1.0)
            end
    end

    @testset "TupleFitnessScheme" begin
            @testset "general interface" begin
                @testset "ParetoFitnessScheme{1}" begin
                        scheme = ParetoFitnessScheme{1}()
                        @test numobjectives(scheme) == 1
                        @test is_minimizing(scheme)
                        @test isnafitness(nafitness(scheme), scheme)
                        @test isequal(nafitness(scheme), (NaN,))
                        @test isnafitness((NaN,), scheme)
                        @test isnafitness((1.0,), scheme) == false
                end

                @testset "ParetoFitnessScheme{1}(is_minimizing=false)" begin
                        scheme = ParetoFitnessScheme{1}(is_minimizing=false)
                        @test numobjectives(scheme) == 1
                        @test is_minimizing(scheme) == false
                        @test isnafitness(nafitness(scheme), scheme)
                        @test isequal(nafitness(scheme), (NaN,))
                        @test isnafitness((NaN,), scheme)
                        @test isnafitness((1.0,), scheme) == false
                end

                @testset "ParetoFitnessScheme{3}()" begin
                        scheme = ParetoFitnessScheme{3}()
                        @test isequal(nafitness(scheme), (NaN,NaN,NaN))
                        @test numobjectives(scheme) == 3
                        @test is_minimizing(scheme)
                        @test isnafitness(nafitness(scheme), scheme)
                        @test isnafitness((NaN,NaN,NaN), scheme)
                        @test isnafitness((1.0,NaN,2.0), scheme)
                        @test isnafitness((1.0,2.0,2.0), scheme) == false
                end
            end

            @testset "hat_compare(..., ParetoFitnessScheme{1}())" begin
                scheme = ParetoFitnessScheme{1}()

                @test_throws MethodError hat_compare((-1.0,), (1.0, 2.0), scheme)

                @test hat_compare((-1.0,), (1.0,), scheme) == -1
                @test hat_compare((0.0,), (0.3,), scheme) == -1
                @test hat_compare((11.3,), (354.65,), scheme) == -1

                @test hat_compare((-1.0,), (-1.0,), scheme) == 0
                @test hat_compare((0.0,), (0.0,), scheme) == 0
                @test hat_compare((0.2,), (0.2,), scheme) == 0

                @test hat_compare((1.0,), (0.0,), scheme) == 1
                @test hat_compare((0.0,), (-0.4,), scheme) == 1
                @test hat_compare((-0.65,), (-34.2,), scheme) == 1
            end

            @testset "hat_compare(..., ParetoFitnessScheme{1}(is_minimizing=false))" begin
                scheme = ParetoFitnessScheme{1}(is_minimizing=false)

                @test hat_compare((-1.0,), (1.0,), scheme) == 1
                @test hat_compare((0.0,), (0.3,), scheme) == 1
                @test hat_compare((11.3,), (354.65,), scheme) == 1

                @test hat_compare((-1.0,), (-1.0,), scheme) == 0
                @test hat_compare((0.0,), (0.0,), scheme) == 0
                @test hat_compare((0.2,), (0.2,), scheme) == 0

                @test hat_compare((1.0,), (0.0,), scheme) == -1
                @test hat_compare((0.0,), (-0.4,), scheme) == -1
                @test hat_compare((-0.65,), (-34.2,), scheme) == -1
            end

            @testset "hat_compare(..., ParetoFitnessScheme{2}(is_minimizing=true))" begin
                scheme = ParetoFitnessScheme{2}(is_minimizing=true)
                @test is_minimizing(scheme)

                @test_throws MethodError hat_compare((-1.0,), (1.0, 2.0), scheme)
                @test_throws MethodError hat_compare((2.0, -1.0, 3.0), (1.0, 2.0), scheme)

                @test hat_compare((-1.0, 0.0), (1.0, 0.0), scheme) == -1
                @test hat_compare((0.0, -1.0), (0.3, 1.0), scheme) == -1

                @test hat_compare((-1.0, 2.0), (-1.0, 2.0), scheme) == 0
                @test hat_compare((0.0, 0.0), (0.0, 0.0), scheme) == 0
                @test hat_compare((0.2, -0.4), (0.2, -0.4), scheme) == 0
                @test hat_compare((0.4, -0.4), (0.2, -0.2), scheme) == 0
                @test hat_compare((-0.65, 1.0), (1.0, -34.2), scheme) == 0

                @test hat_compare((11.3, 100.0), (11.3, 10.0), scheme) == 1
                @test hat_compare((1.0, 0.0), (0.0, -1.0), scheme) == 1
                @test hat_compare((0.0, 34.0), (-0.4, 20.0), scheme) == 1
            end

            @testset "hat_compare(..., ParetoFitnessScheme{2}(is_minimizing=false))" begin
                scheme = ParetoFitnessScheme{2}(is_minimizing=false)
                @test is_minimizing(scheme) == false

                @test hat_compare((-1.0, 0.0), (1.0, 0.0), scheme) == 1
                @test hat_compare((0.0, -1.0), (0.3, 1.0), scheme) == 1

                @test hat_compare((-1.0, 2.0), (-1.0, 2.0), scheme) == 0
                @test hat_compare((0.0, 0.0), (0.0, 0.0), scheme) == 0
                @test hat_compare((0.2, -0.4), (0.2, -0.4), scheme) == 0
                @test hat_compare((0.4, -0.4), (0.2, -0.2), scheme) == 0
                @test hat_compare((-0.65, 1.0), (1.0, -34.2), scheme) == 0

                @test hat_compare((11.3, 100.0), (11.3, 10.0), scheme) == -1
                @test hat_compare((1.0, 0.0), (0.0, -1.0), scheme) == -1
                @test hat_compare((0.0, 34.0), (-0.4, 20.0), scheme) == -1
            end

            @testset "is_better/is_worse/same_fitness(..., ParetoFitnessScheme{2}(is_minimizing=true))" begin
                scheme = ParetoFitnessScheme{2}()

                @test is_better((-1.0, 0.0), (1.0, 0.0), scheme)
                @test is_better((0.0, 0.0), (1.0, 0.0), scheme)
                @test_throws MethodError is_better((0.0,), (1.0,), scheme)

                @test is_better((1.0, 0.0), (-1.0, 0.0), scheme) == false
                @test is_better((1.0, 0.0), (0.0, 0.0), scheme) == false

                @test is_worse((-1.0, 0.0), (1.0, 0.0), scheme) == false
                @test is_worse((0.0, 0.0), (1.0, 0.0), scheme) == false
                @test_throws MethodError is_worse((0.0, 1.0), (1.0,), scheme)

                @test is_worse((1.0, 0.0), (-1.0, 0.0), scheme)
                @test is_worse((1.0, 0.0), (0.0, 0.0), scheme)

                @test same_fitness((-1.0, 0.0), (1.0, 0.0), scheme) == false
                @test same_fitness((0.0, 0.0), (1.0, 0.0), scheme) == false
                @test_throws MethodError same_fitness((0.0,), (1.0,), scheme)

                @test same_fitness((-1.0, 0.0), (-1.0, 0.0), scheme)
                @test same_fitness((0.0, 0.0), (0.0, 0.0), scheme)
                @test_throws MethodError same_fitness((0.0,), (0.0,), scheme)
            end

            @testset "is_better/is_worse/same_fitness(..., ParetoFitnessScheme{2}(is_minimizing=false))" begin
                scheme = ParetoFitnessScheme{2}(is_minimizing=false)

                @test is_better((-1.0, 0.0), (1.0, 0.0), scheme) == false
                @test is_better((0.0, 0.0), (1.0, 0.0), scheme) == false

                @test is_better((1.0, 0.0), (-1.0, 0.0), scheme)
                @test is_better((1.0, 0.0), (0.0, 0.0), scheme)

                @test is_worse((-1.0, 0.0), (1.0, 0.0), scheme)
                @test is_worse((0.0, 0.0), (1.0, 0.0), scheme)

                @test is_worse((1.0, 0.0), (-1.0, 0.0), scheme) == false
                @test is_worse((1.0, 0.0), (0.0, 0.0), scheme) == false

                @test same_fitness((-1.0, 0.0), (1.0, 0.0), scheme) == false
                @test same_fitness((0.0, 0.0), (1.0, 0.0), scheme) == false

                @test same_fitness((-1.0, 0.0), (-1.0, 0.0), scheme)
                @test same_fitness((0.0, 0.0), (0.0, 0.0), scheme)
            end
    end
    @testset "ϵ-box FitnessScheme" begin
            @testset "ϵ-index()" begin
                for (u, ϵ, is_minim, u_ix, delta) in [
                                (2.3, 1.0, true, 2, 0.3),
                                (2.3, 0.1, true, 23, 0.0),
                                (2.3, 1.0, false, 3, 0.7),
                                (1.0, 1.0, true, 1.0, 0.0),
                                (2.0, 2.0, false, 1.0, 0.0),
                                (-1.1, 1.0, true, -2, 0.9),
                                (-5.4, 1.0, false, -5, 0.4)
                        ]
                    res = BlackBoxOptim.ϵ_index(u, ϵ, Val{is_minim})
                    @test res[1] == u_ix
                    @test isapprox(res[2], delta; atol=1E-12)
                end
            end

            @testset "IndexedTupleFitness" begin
                    fit = (0.55, 0.3, -0.21)
                    ifit1 = IndexedTupleFitness(fit, sum(fit), 0.1, Val{true})
                    @test ifit1.orig == fit
                    @test ifit1.agg == sum(fit)
                    @test ifit1.index == (5, 3, -3)
                    @test isapprox(ifit1.dist, norm([0.05, 0.09])/0.1)

                    ifit2 = IndexedTupleFitness(fit, sum(fit), 0.1, Val{false})
                    @test ifit2.orig == fit
                    @test ifit2.agg == sum(fit)
                    @test ifit2.index == (6, 3, -2)
                    @test isapprox(ifit2.dist, norm([0.05, 0.01])/0.1)

                    ifit3 = convert(IndexedTupleFitness, fit,
                                                    EpsBoxDominanceFitnessScheme{3}(0.1))
                    @test ifit3 == ifit1
                    ifit4 = convert(IndexedTupleFitness, fit,
                                                    EpsBoxDominanceFitnessScheme{3}(0.1, is_minimizing=false))
                    @test ifit4 == ifit2
            end

            @testset "hat_compare(..., EpsBoxDominanceFitnessScheme{2}(...))" begin
                minscheme = EpsBoxDominanceFitnessScheme{2}(0.1, is_minimizing=true)
                maxscheme = EpsBoxDominanceFitnessScheme{2}(0.1, is_minimizing=false)

                @test isnafitness(nafitness(minscheme), minscheme)
                @test isnafitness(nafitness(maxscheme), minscheme)
                @test isnafitness(nafitness(minscheme), maxscheme)
                @test isnafitness(nafitness(maxscheme), maxscheme)

                @test hat_compare((-1.0, 0.0), (-1.0, 0.0), minscheme) == (0, true)
                @test hat_compare((-1.0, 0.0), (-1.0, 0.0), maxscheme) == (0, true)
                @test hat_compare((-1.0, 0.0), (1.0, 0.0), minscheme) == (-1, false)
                @test hat_compare((-1.0, 0.0), (1.0, 0.0), maxscheme) == (1, false)
                @test hat_compare((-1.0, 0.0), (-0.99, 0.0), minscheme) == (-1, true)
                @test hat_compare((-1.0, 0.0), (-0.99, 0.0), maxscheme) == (1, false)
                @test hat_compare((-1.0, 0.0), (-1.01, 0.0), minscheme) == (1, false)
                @test hat_compare((-1.0, 0.0), (-1.01, 0.0), maxscheme) == (-1, true)

                @test hat_compare((-0.97, 0.44), (-0.96, 0.43), minscheme) == (0, true)
                @test hat_compare((-0.97, 0.44), (-0.96, 0.43), maxscheme) == (0, true)
                @test hat_compare((-1.0, 0.5), (-1.0, 0.51), minscheme) == (-1, true)
                @test hat_compare((-1.0, 0.5), (-1.0, 0.51), maxscheme) == (1, false)

                @test hat_compare((0.9, 0.0), (1.0, 0.0), minscheme) == (-1, false)
                @test hat_compare((0.9, 0.0), (1.0, 0.0), maxscheme) == (1, false)
                @test hat_compare((0.95, 0.0), (0.98, 0.0), minscheme) == (-1, true)
                @test hat_compare((0.95, 0.0), (0.98, 0.0), maxscheme) == (1, true)

                @test hat_compare((-1.0, 0.05), (1.0, 0.0), minscheme) == (-1, false)
                @test hat_compare((-1.0, 0.05), (1.0, 0.0), maxscheme) == (0, false)
                @test hat_compare((0.9, 0.05), (1.0, 0.0), minscheme) == (-1, false)
                @test hat_compare((0.9, 0.05), (1.0, 0.0), maxscheme) == (0, false)
                @test hat_compare((0.95, 0.05), (0.98, 0.0), minscheme) == (-1, true)
                @test hat_compare((0.95, 0.05), (0.98, 0.0), maxscheme) == (-1, false)
                @test hat_compare((0.95, 0.05), (0.98, 0.01), maxscheme) == (-1, true)
                @test hat_compare((0.95, 0.05), (0.96, 0.0), minscheme) == (1, true)
                @test hat_compare((0.95, 0.05), (0.96, 0.0), maxscheme) == (-1, false)
                @test hat_compare((0.95, 0.05), (0.96, 0.04), minscheme) == (-1, true)
                @test hat_compare((0.95, 0.05), (0.96, 0.04), maxscheme) == (-1, true)

                @test hat_compare((6.9928266604943286,0.03386770153536918),
                                                    (7.007808609410211,0.032833634435236035), minscheme, 0) == (-1, false)
                @test hat_compare((6.9928266604943286,0.03386770153536918),
                                                    (7.007808609410211,0.032833634435236035), minscheme, 1) == (-1, false)
                @test hat_compare((6.9928266604943286,0.03386770153536918),
                                                    (7.007808609410211,0.032833634435236035), minscheme, -1) == (-1, false)
            end
    end
end
