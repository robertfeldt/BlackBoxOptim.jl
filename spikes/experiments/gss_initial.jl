using BlackBoxOptim

require("../generating_set_search.jl")
require("../experiment_framework.jl")

p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], 2^5)
problem = BlackBoxOptim.shifted(p)
ps = {1 => problem}

sf(p) = begin
  generating_set_search(p)
end

@time repeated_runs(sf, ps, 2; experiment = "gss_initial")
