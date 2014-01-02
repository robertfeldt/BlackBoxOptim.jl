using Distributions
using BlackBoxOptim

require("../solis_wets_raw.jl")
require("../experiment_framework.jl")

gaussianrng = (sigma) -> Normal(0, sigma)
levyrng = (sigma) -> Levy(0, sigma)
cauchyrng = (sigma) -> Cauchy(0, sigma)

#p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], 2^3)
#problem = BlackBoxOptim.shifted(p)
#tic()
#xb, fb, nf, r, a = solis_wets(problem; max_fail_steps = 3, 
#  max_evals_per_dim = 1e5, rnggen = cauchyrng, ftol = 1e-7)
#t = toq()

#println("Best solution: $(xb)")
#println("\nTime taken: $(t)")
#println("fevals = $(nf), reason = $(r), fitness = $(fb)")

p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], 2^5)
problem = BlackBoxOptim.shifted(p)
ps = {1 => problem}

sf(p) = begin
  solis_wets(p; max_fail_steps = 3, rnggen = cauchyrng)
end

@time repeated_runs(sf, ps, 2; experiment = "solis_wets_initial")
