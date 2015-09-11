facts("Problems") do

  context("Can create fixed dimensional problem with search space as array") do
    p = BlackBoxOptim.fixeddim_problem((x)->x; search_space = [(0.0, 1.0)])
    @fact typeof(p) --> BlackBoxOptim.FixedDimProblem
  end

  context("Exception if no valid search space supplied when creating fixed dimensional problem") do
    @fact_throws BlackBoxOptim.fixeddim_problem((x)->x; search_space = (0.0, 1.0))
  end

end