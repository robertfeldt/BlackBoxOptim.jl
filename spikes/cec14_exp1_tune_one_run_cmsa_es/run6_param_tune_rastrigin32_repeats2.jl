require("../../src/BlackBoxOptim.jl")
require("cec_cmsa_es.jl")
require("../experiment_framework.jl")
require("../../src/experiments/parameter_experiment.jl")

using BlackBoxOptim

for dim in [32]

  n = dim
  p = BlackBoxOptim.as_fixed_dim_problem(BlackBoxOptim.example_problems["Rastrigin"], n)
  p = BlackBoxOptim.shifted(p)
  diameter = minimum(diameters(search_space(p)))

  pe = ParameterExperiment(
  ["lambda",
   "mu",
   "sigma",
   "decompose_covar_prob",
   "tau",
   "tau_c"
  ],
  [((ds, ps) -> int(10^ds[1])),
   ((ds, ps) -> int(max(1.0, ds[2]*ps[1]))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4]),
   ((ds, ps) -> 1 / sqrt(n * ds[5])),
   ((ds, ps) -> 1 + n * (n+1) / (ds[6]*ps[2]))
  ],
  [( log10(4.0 + log(n)), log10(10*n*n)),
   ( 1/100, 1/2),
   (-3.0, 0.0),
   (0.01, 0.99),
   (0.2, 10.0),
   (0.01, 100.0)
  ]
 )

  for num_repeats in [2]
    prefix = "exp1_6_cmsa_es_6param_log10_$(name(p))$(n)_$(num_repeats)"
    explore_parameters(pe, cmsa_es, p; experiment_prefix = prefix, num_repeats = num_repeats)
  end

end
