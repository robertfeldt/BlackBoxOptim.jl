@testset "Generating set search" begin

    ss = RectSearchSpace(3, (0.0, 1.0))
    @test BlackBoxOptim.calc_initial_step_size(ss) == (0.5 * (1.0 - 0.0))

    ss = RectSearchSpace(3, (-1.2, 42.0))
    @test BlackBoxOptim.calc_initial_step_size(ss, 0.80) == (0.80 * (42.0 + 1.2))  

end
