using Distributions
using BlackBoxOptim

require("../solis_wets.jl")
require("../experiment_framework.jl")

gaussianrng = (sigma) -> Normal(0, sigma)
levyrng = (sigma) -> Levy(0, sigma)
cauchyrng = (sigma) -> Cauchy(0, sigma)

ps = BlackBoxOptim.as_fixed_dim_problem_set(BlackBoxOptim.example_problems, 2^5)

sf(p) = begin
  solis_wets(p; max_fail_steps = 3, rnggen = gaussianrng)
end

@time repeated_runs(sf, ps, 25; experiment = "solis_wets_2")
