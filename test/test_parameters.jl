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

end
