@testset "Single objective functions" begin

@testset "Sphere" begin
    p = BlackBoxOptim.example_problems["Sphere"]
    sphere = objfunc(p)

    @test sphere([0]) == 0
    @test sphere([1]) == 1
    @test sphere([1, 2]) == 5
    @test sphere([1, 2, 3]) == 14
    @test sphere([-1, 2, -3]) == 14
    @test_throws MethodError sphere([])

    p2 = instantiate(p, 3)
    @test numdims(p2) == 3
    @test dimrange(search_space(p2)) == [(-100.0, 100.0), (-100.0, 100.0), (-100.0, 100.0)]
end

@testset "Schwefel2.22" begin
    p = BlackBoxOptim.example_problems["Schwefel2.22"]
    schwefel2_22 = objfunc(p)

    @test schwefel2_22([0]) == 0
    @test schwefel2_22([1]) == 2
    @test schwefel2_22([1, 2]) == (1+2)+(1*2)
    @test schwefel2_22([1, 2, 3]) == (1+2+3)+(1*2*3)
    @test schwefel2_22([-1, 2, -3]) == (1+2+3)+(1*2*3)
    @test_throws MethodError schwefel2_22([])

    p2 = instantiate(p, 4)
    @test numdims(p2) == 4
    @test dimrange(search_space(p2)) == [(-10.0, 10.0), (-10.0, 10.0), (-10.0, 10.0), (-10.0, 10.0)]
end

@testset "Schwefel1.2" begin
    p = BlackBoxOptim.example_problems["Schwefel1.2"]
    schwefel1_2 = objfunc(p)

    @test schwefel1_2([0]) == 0
    @test schwefel1_2([1]) == 1
    @test schwefel1_2([1, 2]) == 1+9
    @test schwefel1_2([1, 2, 3]) == 1+9+36
    @test schwefel1_2([-1, 2, -3]) == 1+1+4
    @test schwefel1_2([]) == 0
end

@testset "Schwefel2.21" begin
    p = BlackBoxOptim.example_problems["Schwefel2.21"]
    schwefel2_21 = objfunc(p)

    @test schwefel2_21([0]) == 0
    @test schwefel2_21([1]) == 1
    @test schwefel2_21([1, 2]) == 2
    @test schwefel2_21([1, 2, 3]) == 3
    @test schwefel2_21([-1, 2, -3]) == 3
    @test_throws MethodError schwefel2_21([])
end

@testset "Rosenbrock" begin
    p = BlackBoxOptim.example_problems["Rosenbrock"]
    rosenbrock = objfunc(p)

    @test rosenbrock([1, 2]) == 100
    @test rosenbrock([1, 2, 3]) == 201
    @test rosenbrock([-1, 2, -3]) == 5005
    @test_throws MethodError rosenbrock([])
end

end
