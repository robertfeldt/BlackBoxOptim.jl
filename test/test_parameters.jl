facts("Parameters") do

  context("When no parameters") do

    ps = Parameters()
    @fact ps[:NotThere] => nothing
    @fact ps["Neither there"] => nothing

    ps[:a] = 1
    @fact ps[:a] => 1

  end

  context("With one parameter in one set") do

    ps = Parameters({:a => 1})

    @fact ps[:a] => 1
    @fact ps["a"] => 1

    @fact ps[:A] => nothing
    @fact ps["A"] => nothing
    @fact ps[:B] => nothing
    @fact ps["B"] => nothing

  end

  context("With parameters in multiple sets") do

    ps = Parameters({:a => 1, "c" => 4}, {:a => 2, :b => 3}, {:c => 5})

    @fact ps[:a] => 1
    @fact ps["a"] => 1

    @fact ps[:c] => 4
    @fact ps["c"] => 4

    @fact ps[:b] => 3
    @fact ps["b"] => 3

    @fact ps[:A] => nothing
    @fact ps["A"] => nothing
    @fact ps[:B] => nothing
    @fact ps["B"] => nothing

  end

  context("Updating parameters after construction") do

    ps = Parameters({:a => 1, "c" => 4}, {:a => 2, :b => 3}, {:c => 5})

    ps[:c] = 6
    ps["b"] = 7

    @fact ps[:a] => 1
    @fact ps["a"] => 1

    @fact ps[:c] => 6
    @fact ps["c"] => 6

    @fact ps[:b] => 7
    @fact ps["b"] => 7

  end

  context("Constructing from another parameters object") do

    ps1 = Parameters({:a => 1, "c" => 4}, {:a => 2, :b => 3})
    ps2 = Parameters({:a => 5}, ps1, {:c => 6})

    @fact ps2[:a] => 5
    @fact ps2[:c] => 4

  end

  context("Get key without default") do

    ps = Parameters({:a => 1, "c" => 4}, {:a => 2, :b => 3})
    @fact get(ps, :a) => 1
    @fact get(ps, :b) => 3
    @fact get(ps, :d) => nothing

  end

  context("Get key without default") do

    ps = Parameters({:a => 1, "c" => 4}, {:a => 2, :b => 3})
    @fact get(ps, :d, 10) => 10

  end

  context("Merge with Parameters or Dict") do

    ps = Parameters({:a => 1, "c" => 4}, {:a => 2, :b => 3})
    ps2 = merge(ps, {:d => 5, :a => 20})
    @fact ps2[:d] => 5
    @fact ps2[:b] => 3
    @fact ps2[:a] => 20

  end
end
