@testset "Random sampling on unit, n-dimensional sphere" begin

    uv1 = BlackBoxOptim.sample_unit_hypersphere(1)
    @test size(uv1, 1) == 1

    uv2 = BlackBoxOptim.sample_unit_hypersphere(2)
    @test size(uv2, 1) == 2

    uv3 = BlackBoxOptim.sample_unit_hypersphere(3, 4)
    @test size(uv3, 1) == 3
    @test size(uv3, 2) == 4

    for i in 1:100
        n = rand(1:100)
        num = rand(2:10)
        u = BlackBoxOptim.sample_unit_hypersphere(n, num)
        @test size(u, 1) == n
        @test size(u, 2) == num
        for j in 1:num
            @test isapprox( norm(u[:,j]), 1.0 )
        end
    end

end

@testset "Random direction generator" begin

    rdg1 = BlackBoxOptim.RandomDirectionGen(2, 3)
    ds1 = BlackBoxOptim.directions_for_k(rdg1, 1)
    @test size(ds1) == (2, 3)

    rdg2 = BlackBoxOptim.RandomDirectionGen(10, 17)
    ds2 = BlackBoxOptim.directions_for_k(rdg2, 1)
    @test size(ds2) == (10, 17)

end

@testset "Mirrored random direction generator" begin

    mrdg1 = BlackBoxOptim.MirroredRandomDirectionGen(2, 4)
    ds1 = BlackBoxOptim.directions_for_k(mrdg1, 1)
    @test size(ds1) == (2, 4)
    @test ds1[:,3] == -ds1[:,1]
    @test ds1[:,4] == -ds1[:,2]

    mrdg2 = BlackBoxOptim.MirroredRandomDirectionGen(10, 6)
    ds2 = BlackBoxOptim.directions_for_k(mrdg2, 1)
    @test size(ds2) == (10, 6)
    @test ds2[:,4] == -ds2[:,1]
    @test ds2[:,5] == -ds2[:,2]
    @test ds2[:,6] == -ds2[:,3]

    # Must be even number of directions
    @test_throws ArgumentError BlackBoxOptim.MirroredRandomDirectionGen(10, 1)
    @test_throws ArgumentError BlackBoxOptim.MirroredRandomDirectionGen(10, 3)
    @test_throws ArgumentError BlackBoxOptim.MirroredRandomDirectionGen(10, 7)

end
