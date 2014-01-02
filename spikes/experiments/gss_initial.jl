using BlackBoxOptim

require("../generating_set_search.jl")
require("../experiment_framework.jl")

sf(p) = begin
  generating_set_search(p)
end

ps = BlackBoxOptim.as_fixed_dim_problem_set(BlackBoxOptim.example_problems, 2.^[2,6])

@time repeated_runs(sf, ps, 3; experiment = "gss")
