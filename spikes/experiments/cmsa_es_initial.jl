using BlackBoxOptim

require("../covar_matrix_self_adaptation_es.jl")
require("../experiment_framework.jl")

ps = BlackBoxOptim.as_fixed_dim_problem_set(BlackBoxOptim.example_problems, 2^5)

sf(p) = begin
  p = BlackBoxOptim.shifted(p)
  n = numdims(p)
  lambda = 4 + 3 * int(ceil(log10(n)))
  mu = 2 # int(ceil(lambda/2))
  cmsa_es(p, mu = mu, lambda = lambda)
end

@time repeated_runs(sf, ps, 25; experiment = "cmsa_es_2")
