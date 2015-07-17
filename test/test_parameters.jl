facts("DictChain") do
  context("Matching keys and type parameters") do
    dc1 = DictChain{Int,ASCIIString}()

    @fact_throws (dc1[:NotThere] = 3) KeyError
    @fact_throws (dc1[:NotThere] = "abc") KeyError
    @fact_throws (dc1[3] = 3) MethodError
    dc1[3] = "abc"
    @fact dc1[3] => "abc"
  end

  context("merging and chaining") do
    dc1 = DictChain{Symbol,ASCIIString}()

    # incompatible dictionary types
    context("incompatible dictionary types") do
      dc2 = DictChain{Int,ASCIIString}()
      @fact_throws chain(dc1, dc2) MethodError
      @fact_throws merge(dc1, dc2) MethodError
      @fact_throws merge!(dc1, dc2) MethodError

      dc3 = DictChain{Symbol,Int}()
      @fact_throws chain(dc1, dc3) MethodError
      @fact_throws merge(dc1, dc3) MethodError
      @fact_throws merge!(dc1, dc3) MethodError
    end

    d1 = @compat Dict{Symbol,Int}(:a => 1)
    d2 = @compat Dict{Symbol,Int}(:a => 2, :b => 4)
    d3 = @compat Dict{Symbol,Int}(:b => 3, :c => 5)

    context("using constructor") do
      dc = DictChain(d1, d2, d3)
      @fact typeof(dc) => DictChain{Symbol,Int}
      @fact dc[:a], dc[:b], dc[:c] => 1, 4, 5

      dc = DictChain(d2, d1, d3)
      @fact dc[:a], dc[:b], dc[:c] => 2, 4, 5

      dc = DictChain(DictChain(d3, d1), d2)
      @fact dc[:a], dc[:b], dc[:c] => 1, 3, 5
    end

    context("using merge()") do
      dc = merge(merge(d1, d2), d3)
      @fact dc[:a], dc[:b], dc[:c] => 2, 3, 5

      dc = merge(d2, merge(d1, d3))
      @fact dc[:a], dc[:b], dc[:c] => 1, 3, 5

      dc = merge(merge(d3, d1), d2)
      @fact dc[:a], dc[:b], dc[:c] => 2, 4, 5
    end

    context("using merge!()") do
      dc = DictChain{Int,ASCIIString}()
      @fact_throws merge!(dc, d1) MethodError # incompatible types

      dc = DictChain{Symbol,Int}()
      merge!(dc, d1)
      @fact dc[:a] => 1
      @fact_throws dc[:b] KeyError
      merge!(dc, d2)
      @fact dc[:a], dc[:b] => 2, 4
      merge!(dc, d3)
      @fact dc[:a], dc[:b], dc[:c] => 2, 3, 5
    end

    context("using chain()") do
      dc = chain(chain(d1, d2), d3)
      @fact typeof(dc) => DictChain{Symbol,Int}
      @fact dc[:a], dc[:b], dc[:c] => 2, 3, 5

      dc = chain(d2, chain(d1, d3))
      @fact dc[:a], dc[:b], dc[:c] => 1, 3, 5

      dc = chain(chain(d3, d1), d2)
      @fact dc[:a], dc[:b], dc[:c] => 2, 4, 5
    end
  end

  context("converting to Dict") do
    d1 = @compat Dict{Symbol,Int}(:a => 1)
    d2 = @compat Dict{Symbol,Int}(:a => 2, :b => 4)
    d3 = @compat Dict{Symbol,Int}(:a => 3, :b => 5)

    dc = DictChain(d1, d2, d3)
    d123 = convert(Dict{Symbol,Int}, dc)
    @fact typeof(d123) => Dict{Symbol,Int}
    @fact length(d123) => 2
    @fact d123[:a] => 1
    @fact d123[:b] => 4
  end
end

facts("Parameters") do

  context("When no parameters or key type doesn't match") do

    ps = Parameters()
    @fact_throws ps[:NotThere] KeyError
    @fact_throws ps["Neither there"] MethodError

    ps[:a] = 1
    @fact ps[:a] => 1

  end

  context("With one parameter in one set") do

    ps = Parameters(@compat Dict{Symbol,Any}(:a => 1))

    @fact ps[:a] => 1
    @fact_throws ps["a"] MethodError # incorrect key type

    @fact_throws ps[:A] KeyError
    @fact_throws ps[:B] KeyError

  end

  context("With parameters in multiple sets") do

    ps = Parameters(@compat(Dict{Symbol,Any}(:a => 1, :c => 4)),
                    @compat(Dict{Symbol,Any}(:a => 2, :b => 3)),
                    @compat(Dict{Symbol,Any}(:c => 5)))

    @fact ps[:a] => 1

    @fact ps[:c] => 4

    @fact ps[:b] => 3

    @fact_throws ps[:A] KeyError
    @fact_throws ps[:B] KeyError
  end

  context("Updating parameters after construction") do

    ps = Parameters(@compat(Dict{Symbol,Any}(:a => 1, :c => 4)),
                    @compat(Dict{Symbol,Any}(:a => 2, :b => 3)),
                    @compat(Dict{Symbol,Any}(:c => 5)))

    ps[:c] = 6
    ps[:b] = 7

    @fact ps[:a] => 1

    @fact ps[:c] => 6

    @fact ps[:b] => 7

  end

  context("Constructing from another parameters object") do

    ps1 = Parameters(@compat(Dict{Symbol,Any}(:a => 1, :c => 4)),
                     @compat(Dict{Symbol,Any}(:a => 2, :b => 3)))
    ps2 = Parameters(@compat(Dict{Symbol,Any}(:a => 5)), ps1,
                     @compat(Dict{Symbol,Any}(:c => 6)))

    @fact ps1[:a] => 1
    @fact ps2[:a] => 5
    @fact ps2[:c] => 4
  end

  context("Get key without default") do

    ps = Parameters(@compat(Dict{Symbol,Any}(:a => 1, :c => 4)),
                    @compat(Dict{Symbol,Any}(:a => 2, :b => 3)))
    @fact get(ps, :a) => 1
    @fact get(ps, :b) => 3
    @fact get(ps, :d) => nothing

  end

  context("Get key with default") do

    ps = Parameters(@compat(Dict{Symbol,Any}(:a => 1, :c => 4)),
                    @compat(Dict{Symbol,Any}(:a => 2, :b => 3)))
    @fact get(ps, :d, 10) => 10

  end

  context("Merge with Parameters or Dict") do

    ps = Parameters(@compat(Dict{Symbol,Any}(:a => 1, :c => 4)),
                    @compat(Dict{Symbol,Any}(:a => 2, :b => 3)))
    ps2 = chain(ps, @compat(Dict{Symbol,Any}(:d => 5, :a => 20)))
    @fact ps2[:d] => 5
    @fact ps2[:b] => 3
    @fact ps2[:a] => 20

  end
end
