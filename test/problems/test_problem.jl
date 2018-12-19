fsabs(x) = sum(abs, x)
fsum_abs_and_sq(x) = (sum(abs, x), sum(abs2, x))

@testset "minimization_problem()" begin

@testset "1-dimensional, single-objective sum(abs, x)" begin
    p = minimization_problem(fsabs, "sumabs", (-1.0, 1.0), 1, 0.0)

    @test fitness_scheme(p) == MinimizingFitnessScheme
    @test numobjectives(p) == 1
    @test numdims(p) == 1
    @test search_space(p) == RectSearchSpace(1, (-1.0, 1.0))

    @test fitness([0.0], p) == 0.0
    @test fitness([1.2], p) == 1.2
    @test fitness([-1.9], p) == 1.9
end

@testset "3-dimensional, single-objective sum(abs, x)" begin
    p = minimization_problem(fsabs, "sumabs", (-1.0, 1.0), 3, 0.0)

    @test fitness_scheme(p) == MinimizingFitnessScheme
    @test numobjectives(p) == 1
    @test numdims(p) == 3
    @test search_space(p) == RectSearchSpace(3, (-1.0, 1.0))

    @test fitness([0.0, 1.0, 2.0], p) == 3.0
    @test fitness([-1.0, 1.0, 2.0], p) == 4.0
end

@testset "1-dimensional, multi-objective sumabs_sumsq" begin
    ss = RectSearchSpace(1)
    p = FunctionBasedProblem(fsum_abs_and_sq, "sumabs_sumsq",
                             ParetoFitnessScheme{2}(), ss, (0.0, 0.0))

    @test numdims(p) == 1
    @test numobjectives(p) == 2
    @test search_space(p) == ss

    @test fitness([0.0], p) == (0.0, 0.0)
    @test fitness([1.2], p) == (1.2, 1.2^2)
    @test fitness([-1.9], p) == (1.9, 1.9^2)
end

@testset "MinimizationProblemFamily" begin
    pfam = MinimizationProblemFamily(fsabs, "sumabs", (0.0, 1.0), 0.0)
    p1 = instantiate(pfam, 1)

    @test numobjectives(p1) == 1
    @test numdims(p1) == 1

    @test fitness([0.0], p1) == 0.0
    @test fitness([1.2], p1) == 1.2
    @test fitness([-1.9], p1) == 1.9

    p3 = instantiate(pfam, 3)

    @test numobjectives(p3) == 1
    @test numdims(p3) == 3

    @test fitness([0.0, 1.0, 2.0], p3) == 3.0
    @test fitness([-1.0, 1.0, 2.0], p3) == 4.0
end

end

@testset "ShiftedAndBiasedProblem" begin

@testset "Only x shifted 1-dim" begin
    p = minimization_problem(fsabs, "sumabs", (-1.0, 1.0), 1, 0.0)
    sp = BlackBoxOptim.ShiftedAndBiasedProblem(p, xshift=[0.5])

    @test BlackBoxOptim.orig_problem(sp) === p
    @test numobjectives(sp) == 1
    @test numdims(sp) == 1
    @test search_space(sp) == RectSearchSpace(1, (-1.0, 1.0))
    xs = sp.xshift
    @test xs == [0.5]
    @test sp.fitshift == 0.0

    @test fitness([0.0], sp) == 0.5
    @test fitness([1.2], sp) == 0.7
    @test fitness([-1.9], sp) == 2.4
end

@testset "Shifted and biased 2-dim" begin
    ss = RectSearchSpace(2, (-0.5, 1.0))
    p = minimization_problem(fsabs, "sumabs", (-0.5, 1.0), 2, 0.0)
    sp = BlackBoxOptim.ShiftedAndBiasedProblem(p; fitshift = 1.3)

    @test BlackBoxOptim.orig_problem(sp) === p
    @test numobjectives(sp) == 1
    @test numdims(sp) == 2
    @test sp.xshift == [0.0, 0.0]
    @test sp.fitshift == 1.3
    @test search_space(sp) == RectSearchSpace(2, (-0.5, 1.0))

    xs = sp.xshift
    @test fitness([0.0, 1.0], sp) == 1.0 + 1.3
    @test fitness([1.2, -1.3], sp) == 2.5 + 1.3
    @test fitness([-1.9, 0.0], sp) == 1.9 + 1.3
end

@testset "Within ftol" begin
    subp = minimization_problem(fsabs, "sumabs", (-1.0, 1.0), 2, 0.0)

    @test fitness_is_within_ftol(subp, 0.1, 0.2)
    @test !fitness_is_within_ftol(subp, 0.1, 0.09)
end

end
