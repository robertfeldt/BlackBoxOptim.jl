facts("Assign ranks within tolerance") do

  context("Ranks correctly if none are within tolerance of each other") do
    res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([1.0, 2.0, 3.0])
    @fact length(res) => 3
    @fact res[1] => (1, 1.0, 1.0)
    @fact res[2] => (2, 2.0, 2.0)
    @fact res[3] => (3, 3.0, 3.0)

    res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([3.0])
    @fact length(res) => 1
    @fact res[1] => (1, 3.0, 3.0)
  end

  context("Ranks correctly if some are within tolerance of each other") do
    res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([1.0, 1.0+1e-6, 200.0]; tolerance = 1e-5)
    @fact length(res) => 3
    @fact res[1] => (1, 1.0, 1.0)
    @fact res[2] => (1, 1.0+1e-6, 1.0+1e-6)
    @fact res[3] => (3, 200.0, 200.0)
  end

  context("Ranks correctly if all are within tolerance of each other") do
    res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([1.0, 1.0+1e-6, 1.0-1e-6]; tolerance = 1e-5)
    @fact length(res) => 3
    @fact res[1] => (1, 1.0-1e-6, 1.0-1e-6)
    @fact res[2] => (1, 1.0, 1.0)
    @fact res[3] => (1, 1.0+1e-6, 1.0+1e-6)
  end

  context("Ranks in reverse if none are within tolerance of each other") do
    res = BlackBoxOptim.Utils.assign_ranks_within_tolerance([1.0, 20.0, 13.0]; rev = true)
    @fact length(res) => 3
    @fact res[1] => (1, 20.0, 20.0)
    @fact res[2] => (2, 13.0, 13.0)
    @fact res[3] => (3, 1.0, 1.0)
  end

  context("Ranks correctly for a complex example testing many aspects") do
    values = shuffle([(-11.0, :a), (1.0, :b), (1.0+1e-6, :c), (1.0-1e-6, :d), 
      (2000, :e), (3000, :f)])
    res = BlackBoxOptim.Utils.assign_ranks_within_tolerance(values; 
      tolerance = 1e-5, by = (p) -> p[1], rev = true)
    @fact length(res) => length(values)
    @fact res[1] => (1, (3000, :f), 3000)
    @fact res[2] => (2, (2000, :e), 2000)
    @fact res[3] => (3, (1.0+1e-6, :c), 1.0+1e-6)
    @fact res[4] => (3, (1.0, :b), 1.0)
    @fact res[5] => (3, (1.0-1e-6, :d), 1.0-1e-6)
    @fact res[6] => (6, (-11.0, :a), -11.0)
  end

end