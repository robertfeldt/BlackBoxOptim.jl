@testset "DictChain" begin
    @testset "Matching keys and value types for get()/set() methods" begin
        dc1 = DictChain{Int,String}()

        @test_throws KeyError dc1[:NotThere]
        @test_throws MethodError (dc1[:NotThere] = 3)
        @test_throws MethodError (dc1[:NotThere] = "abc")
        @test_throws MethodError (dc1[3] = 3)

        dc1[3] = "abc"
        @test dc1[3] == "abc"
    end

    @testset "merging and chaining" begin
        dc1 = DictChain{Symbol,String}()

        # incompatible dictionary types
        @testset "incompatible dictionary types" begin
            dc2 = DictChain{Int,String}()
            @test_throws MethodError chain(dc1, dc2)
            @test_throws MethodError merge(dc1, dc2)
            @test_throws MethodError merge!(dc1, dc2)

            dc3 = DictChain{Symbol,Int}()
            @test_throws MethodError chain(dc1, dc3)
            @test_throws MethodError merge(dc1, dc3)
            @test_throws MethodError merge!(dc1, dc3)
        end

        d1 = Dict{Symbol,Int}(:a => 1)
        d2 = Dict{Symbol,Int}(:a => 2, :b => 4)
        d3 = Dict{Symbol,Int}(:b => 3, :c => 5)

        @testset "using constructor" begin
            dc = DictChain(d1, d2, d3)
            @test typeof(dc) == DictChain{Symbol,Int}
            @test (dc[:a], dc[:b], dc[:c]) == (1, 4, 5)

            dc = DictChain(d2, d1, d3)
            @test (dc[:a], dc[:b], dc[:c]) == (2, 4, 5)

            dc = DictChain(DictChain(d3, d1), d2)
            @test (dc[:a], dc[:b], dc[:c]) == (1, 3, 5)
        end

        @testset "using merge()" begin
            dc = merge(DictChain(d2, d1), d3)
            @test (dc[:a], dc[:b], dc[:c]) == (2, 3, 5)

            dc = merge(d2, DictChain(d3, d1))
            @test (dc[:a], dc[:b], dc[:c]) == (1, 3, 5)

            dc = merge(DictChain(d1, d3), d2)
            @test (dc[:a], dc[:b], dc[:c]) == (2, 4, 5)
        end

        @testset "using merge!()" begin
            dc = DictChain{Int,String}()
            @test_throws MethodError merge!(dc, d1) # incompatible

            dc = DictChain{Symbol,Int}()
            merge!(dc, d1)
            @test dc[:a] == 1
            @test_throws KeyError dc[:b]
            merge!(dc, d2)
            @test (dc[:a], dc[:b]) == (2, 4)
            merge!(dc, d3)
            @test (dc[:a], dc[:b], dc[:c]) == (2, 3, 5)
        end

        @testset "using chain()" begin
            dc = chain(chain(d1, d2), d3)
            @test typeof(dc) == DictChain{Symbol,Int}
            @test (dc[:a], dc[:b], dc[:c]) == (2, 3, 5)

            dc = chain(d2, chain(d1, d3))
            @test (dc[:a], dc[:b], dc[:c]) == (1, 3, 5)

            dc = chain(chain(d3, d1), d2)
            @test (dc[:a], dc[:b], dc[:c]) == (2, 4, 5)

            dc = chain(d1, d2, d3)
            @test (dc[:a], dc[:b], dc[:c]) == (2, 3, 5)

            dc = chain(d3, d1, d2)
            @test (dc[:a], dc[:b], dc[:c]) == (2, 4, 5)
        end

        @testset "Check chaining from real use case" begin
            # This sequence is used in compare_optimizers.jl:repeated_bboptimize
            parameters = BlackBoxOptim.EMPTY_PARAMS
            ftol = 1e-5
            params = chain(parameters, ParamsDict(:FitnessTolerance => ftol))
            @test params[:FitnessTolerance] == ftol
        end
    end

    @testset "converting to Dict" begin
        d1 = Dict{Symbol,Int}(:a => 1)
        d2 = Dict{Symbol,Int}(:a => 2, :b => 4)
        d3 = Dict{Symbol,Int}(:a => 3, :b => 5)

        dc = DictChain(d1, d2, d3)
        d123 = convert(Dict{Symbol,Int}, dc)
        @test typeof(d123) == Dict{Symbol,Int}
        @test length(d123) == 2
        @test d123[:a] == 1
        @test d123[:b] == 4
    end

    @testset "show()" begin
        d1 = Dict{Symbol,Int}(:a => 1)
        d2 = Dict{Symbol,Int}(:a => 2, :b => 4)
        d3 = Dict{Symbol,Int}(:a => 3, :b => 5)

        dc = DictChain(d1, d2, d3)
        iob = IOBuffer()
        show(iob, dc)
        @test replace(String(take!(iob)), ' '=>"") == "BlackBoxOptim.DictChain{Symbol,$Int}[Dict(:a=>1),Dict(:a=>2,:b=>4),Dict(:a=>3,:b=>5)]"
    end

    @testset "flatten" begin
        d1 = Dict{Symbol,Int}(:a => 1)
        d2 = Dict{Symbol,Int}(:a => 2, :b => 4)
        d3 = Dict{Symbol,Int}(:a => 3, :b => 5)

        dc = DictChain(d1, d2, d3)

        fd = flatten(dc)
        @test fd[:a] == 1
        @test fd[:b] == 4
        @test sort(collect(keys(fd))) == [:a, :b]
    end
end

@testset "Parameters" begin

    @testset "When no parameters or key type doesn't match" begin
        ps = ParamsDictChain()
        @test isa(ps, Parameters)
        @test_throws KeyError ps[:NotThere]
        @test_throws KeyError ps["Neither there"]

        ps[:a] = 1
        @test ps[:a] == 1
    end

    @testset "With one parameter in one set" begin
        ps = ParamsDictChain(ParamsDict(:a => 1))
        @test isa(ps, Parameters)

        @test ps[:a] == 1
        @test_throws KeyError ps["a"] # incorrect key

        @test_throws KeyError ps[:A]
        @test_throws KeyError ps[:B]
    end

    @testset "With parameters in multiple sets" begin
        ps = ParamsDictChain(ParamsDict(:a => 1, :c => 4),
                             ParamsDict(:a => 2, :b => 3),
                             ParamsDict(:c => 5))
        @test isa(ps, Parameters)

        @test ps[:a] == 1
        @test ps[:c] == 4
        @test ps[:b] == 3

        @test_throws KeyError ps[:A]
        @test_throws KeyError ps[:B]
    end

    @testset "Updating parameters after construction" begin
        ps = ParamsDictChain(ParamsDict(:a => 1, :c => 4),
                             ParamsDict(:a => 2, :b => 3),
                             ParamsDict(:c => 5))

        ps[:c] = 6
        ps[:b] = 7

        @test ps[:a] == 1
        @test ps[:c] == 6
        @test ps[:b] == 7
    end

    @testset "Constructing from another parameters object" begin
        ps1 = ParamsDictChain(ParamsDict(:a => 1, :c => 4),
                              ParamsDict(:a => 2, :b => 3))
        ps2 = ParamsDictChain(ParamsDict(:a => 5), ps1,
                              ParamsDict(:c => 6))

        @test ps1[:a] == 1
        @test ps2[:a] == 5
        @test ps2[:c] == 4
    end

    @testset "Get key without default" begin
        ps = ParamsDictChain(ParamsDict(:a => 1, :c => 4),
                             ParamsDict(:a => 2, :b => 3))
        @test get(ps, :a) == 1
        @test get(ps, :b) == 3
        @test get(ps, :d) == nothing
    end

    @testset "Get key with default" begin
        ps = ParamsDictChain(ParamsDict(:a => 1, :c => 4),
                             ParamsDict(:a => 2, :b => 3))
        @test get(ps, :d, 10) == 10
    end

    @testset "Merge with Parameters or Dict" begin
        ps = ParamsDictChain(ParamsDict(:a => 1, :c => 4),
                             ParamsDict(:a => 2, :b => 3))
        ps2 = chain(ps, ParamsDict(:d => 5, :a => 20))
        @test ps2[:d] == 5
        @test ps2[:b] == 3
        @test ps2[:a] == 20
    end
end
