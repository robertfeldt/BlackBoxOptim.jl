require("../../src/BlackBoxOptim.jl")
require("../covar_matrix_self_adaptation_es.jl")
require("../experiment_framework.jl")
require("../../src/experiments/parameter_experiment.jl")

using BlackBoxOptim

for dim in [2,4,8,16,32,64]

  n = dim
  p = BlackBoxOptim.as_fixed_dim_problem(BlackBoxOptim.example_problems["Sphere"], n)
  p = BlackBoxOptim.shifted(p)
  diameter = minimum(diameters(search_space(p)))

  pe = ParameterExperiment(
  ["lambda", "mu", "sigma", "decompose_covar_prob",
   "max_successes_before_increasing", "max_failures_before_decreasing",
   "max_rounds_without_improvement", "tau_c"
  ],
  [((ds, ps) -> int(10^ds[1])),
   ((ds, ps) -> int(max(1.0, ds[2]*ps[1]))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4]),
   ((ds, ps) -> int(ds[5])),
   ((ds, ps) -> int(ds[6])),
   ((ds, ps) -> int(10^ds[7])),
   ((ds, ps) -> 1 + n * (n+1) / (ds[8]*ps[2]))
  ],
  [( log10(4.0 + log(n)), log10(4*n*n)),
   ( 1/50, 1/2),
   (-3.0, 0.0),
   (0.01, 0.99),
   (-1000.0, 1000.0), # Inactive if less than 1
   (-1000.0, 1000.0), # Inactive if less than 1
   (1.0, 4.0),
   (0.01, 100.0)
  ]
 )

 # Do when with fewer repeats just to test things...
 explore_parameters(pe, cmsa_es, p; experiment_prefix = "cmsa_es", num_repeats = 10)

end
