using BlackBoxOptim

require("../covar_matrix_self_adaptation_es.jl")
require("../experiment_framework.jl")
require("../../src/experiments/parameter_experiment.jl")

n = 8
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe = cmsa_es_exp2 = ParameterExperiment(
  ["covar_learning_rate", "sigma", "lambda", "mu",
    "covarMatrixSampler", "utilitiesFunc"
  ],
  [((ds, ps) -> ds[1]),
   ((ds, ps) -> diameter * 10^ds[2]),
   ((ds, ps) -> int(round(ds[3]))),
   ((ds, ps) -> int(max(1.0, ceil(ps[3] / ds[4])))),
   ((ds, ps) -> ds[5] < 0.5 ? EigenCovarSampler : CholeskyCovarSampler),
   ((ds, ps) -> ds[6] < 0.5 ? linear_utilities : log_utilities)
  ],
  [( 0.6, 1.0),
   (-2.0, 0.0),
   ( 8.0, float(n*n)),
   ( 2.0, 8.0),
   ( 0.0, 1.0),
   ( 0.0, 1.0)
  ]
)

# It seems clear that (in order of sensitivity)
#  1. covar_learning_rate needs to be high (highest sensitivity)
#  2. mu should be a rather high divisor, at least 4, maybe even 8
#  3. we should use log_utilities (but large variation on this one so unclear)
#  4. sigma should be rather low (but low sensitivity)
#  5. lambda should be high, rather than low (low s)
#  6. the covar sampler does not have a large effect but Cholesky seems somewhat better
#
# We design a new experiment in line with this:

pe3 = cmsa_es_exp3 = ParameterExperiment(
  ["covar_learning_rate", "sigma", "lambda", "mu",
  ],
  [((ds, ps) -> ds[1]),
   ((ds, ps) -> diameter * 10^ds[2]),
   ((ds, ps) -> int(round(ds[3]))),
   ((ds, ps) -> int(max(1.0, ceil(ps[3] / ds[4]))))
  ],
  [( 0.75, 1.0),
   (-4.0, -1.0),
   (float(4*n), float(4*n*n)),
   ( 3.0, 16.0)
  ]
)

# There is quite a lot of variation for all of these covariates, except for
# low values of covariate 1 but conclusions are:
#   1. mu divisor should be larger than 5.6 (3+(16-3)*0.20) and higher values seem to be good
#   2. lambda should be between 0.10-0.35 of the range 4*n to 4*n*n
#   3. covar_learning_rate should be between 0.35-0.60 of the range 0.75-1.0 so in range 0.84-0.90
#   4. sigma should be 0.80-1.0 in its range so 10^[-2, -1] of diameter with higher values giving quicker responses.

# Based on these 2 experiments (exp2 and exp3) we conclude that:
#   1. covar_learning_rate should be in 0.85-0.90
#   2. lambda should be circa 4*n+.20*(4*n*n - 4*n)
#   3. mu divisor should be in 5-16
#   4. sigma should be 0.01-0.10 of the diameter
# Lets now try these settings on the Rosenbrock with 16 dims.

n = 16
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe4 = cmsa_es_exp4 = ParameterExperiment(
  ["covar_learning_rate", "sigma", "lambda", "mu",
  ],
  [((ds, ps) -> ds[1]),
   ((ds, ps) -> diameter * 10^ds[2]),
   ((ds, ps) -> int(round(ds[3]))),
   ((ds, ps) -> int(max(1.0, ceil(ps[3] / ds[4]))))
  ],
  [( 0.84, 0.99),
   (-2.0, -0.5),
   (float(4*n), float(4*n*n)),
   ( 5.0, 20.0)
  ]
)

outfile = "cmsa_es_exp4.csv"
write_csv_header_to_file(pe4, outfile)

# Optimize to minimize the fail rate, i.e. inverted success rate:
# /usr/bin/Rscript ../../R/parameter_experiment.R 4 5 -11 4 cmsa_es_exp4.csv new_runs.json

run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(pe4, outfile, 5; num_repeats = 5)

# /usr/bin/Rscript ../../R/parameter_experiment.R 4 5 -11 4 cmsa_es_exp4.csv new_runs.json

# For Rosenbrock 16 we saw good performance for the same param values as for
# the 8-d version. In detail we conclude:
#   1. best performance for lambda closer to 4*n than higher (highest sens)
#   2. best performance for low mu all the way down to divisor 20 (high sens)
#   3. covar learning rate should be from around 0.89-0.95 but worse if higher
#   4. fairly insensitive to sigma as long as its around 0.05-0.20 of diameter (lowest sens)
#
# Lets now try these settings on the Rosenbrock with 32 dims.

n = 32
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe5 = ParameterExperiment(
  ["covar_learning_rate", "sigma", "lambda", "mu",
  ],
  [((ds, ps) -> ds[1]),
   ((ds, ps) -> diameter * 10^ds[2]),
   ((ds, ps) -> int(round(ds[3]))),
   ((ds, ps) -> int(max(1.0, ceil(ps[3] / ds[4]))))
  ],
  [( 0.84, 0.99),
   (-2.0, -0.5),
   (float(2*n), float(2*n*n)),
   ( 5.0, 25.0)
  ]
)

outfile = "cmsa_es_exp5.csv"
write_csv_header_to_file(pe5, outfile)

# /usr/bin/Rscript ../../R/parameter_experiment.R 4 5 -11 4 cmsa_es_exp5.csv new_runs.json

run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(pe5, outfile, 5; num_repeats = 5)

# Meta-strategy: Not good to run 5 initial runs and then do 25 ones based on EI =>
# only 10% or so for n=32 had high success rate. This is in line with common sense, i.e.
# when doing adaptive sampling one should do it in very small batches (1 run at a time)
# to learn as much as possible about the response before each run).

# For Rosenbrock 32 we saw good performance but for somewhat different values than
# for 8 and 16:
# After 35 runs:
# best: {:covar_learning_rate=>0.9851933684393066,:lambda=>331,:sigma=>0.937811910301762,:mu=>15}
# top5_max: {:covar_learning_rate=>0.9860690669022734,:lambda=>501,:sigma=>2.650768642586489,:mu=>22}
# top5_mean: {:covar_learning_rate=>0.9822724531856971,:lambda=>344,:sigma=>1.4712605685842668,:mu=>18}
# top5_min: {:covar_learning_rate=>0.9780927470649825,:lambda=>168,:sigma=>0.8147138415085076,:mu=>12}
#
#   1. best performance for lambda around 5*n to 15*n but better with lower (highest sens)
#   2. covar learning rate should be from around 0.97-0.98 (high sens)
#   3. best performance for low mu all the way down to divisor 20 (lower sens)
#   4. fairly insensitive to sigma as long as its around 0.02-0.07 of diameter (lowest sens)

# Before addressing Rosenbrock 64-D we added stopping if no progress in 1000*lambda
# func evals (for cmsa_es). I also fixed a bug in the cmsa_es code so that the sigma
# does not increase indefinetely to compensate for bad scaling of the covar matrix update.
# Because of this we also changed back to the tau equations of Beyer et al.
# And added a probability for decomposing the covar, to save time. But it might affect quality.
# So I guess I might need to have to rerun all experiments... ;)

n = 64
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe6 = ParameterExperiment(
  ["tau", "lambda", "mu", "tau_c", "sigma", "decompose_covar_prob"
  ],
  [((ds, ps) -> (2*n)^ds[1]),
   ((ds, ps) -> n*int(round(ds[2]))),
   ((ds, ps) -> int(max(1.0, ceil(ps[2] / ds[3])))),
   ((ds, ps) -> 1 + n * (n+1) / (ds[4]*ps[3])),
   ((ds, ps) -> diameter * 10^ds[5]),
   ((ds, ps) -> ds[6])
  ],
  [(-2.0, 0.0),
   ( 1.0, float(n)),
   ( 4.0, 20.0),
   (0.01, 3.0),
   (-3.0, 0.0),
   (0.01, 1.0)
  ]
)

# First we minimize the magnitude class with ei => 14
outfile = "cmsa_es_exp6.csv"
write_csv_header_to_file(pe6, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 6 7 14 6 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe6, outfile; num_repeats = 1)
for i in 1:33
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 6 1 14 6 $(outfile) new_runs.json sa min`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe6, outfile; num_repeats = 2)
end
# Then we maximize the success rate and select the quickest among the top 10.
for i in 1:10
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 6 1 14 6 $(outfile) new_runs.json sa minquick`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe6, outfile; num_repeats = 2)
end

# Conclusions for finding good solutions (in sensitivity order):
#  1. lambda should be a fairly low multiple of n, 1-8
#  2. values of tau seem influential but lots of variation in its effect so hard to give guidance => fix to theoretic value!?
#  3. tau_c fairly influential and should probably be a divisor of 3*mu rather than 2*mu => fix to 2*mu?
#  4. mu somewhat influential and should rather be a high than a low divisor, maybe in range 8-20 rather than 4-8 of lambda.
#  5. decompose_covar_prob not so influential on its own but interacts with the others higher values seems generally somewhat better, but not strongly and affects time behavior.
#  6. sigma not so important but rather higher (1*diameter) than lower (0.001*diameter)

# Lets do a final Rosenbrock 128-D before going on to other problems.

n = 128
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe7 = ParameterExperiment(
  ["lambda", "mu", "sigma", "decompose_covar_prob"
  ],
  [((ds, ps) -> int(round(n*ds[1]))),
   ((ds, ps) -> int(max(1.0, ceil(ps[1] / ds[2])))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4])
  ],
  [( 0.5, float(n)),
   ( 4.0, 20.0),
   (-2.0, 0.0),
   (0.01, 1.0)
  ]
)

outfile = "cmsa_es_exp7.csv"
write_csv_header_to_file(pe7, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 4 5 10 4 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe7, outfile; num_repeats = 1)

# It seems that covergence is slow. I posit that this may be because of mu being too high =>
# increase size of max divisor in experiments.
# We also try to reduce the decompose_covar_prob max value since the cost of
# decomposing the covar matrix is high for high-dimensional problems.

pe8 = ParameterExperiment(
  ["lambda", "mu", "sigma", "decompose_covar_prob"
  ],
  [((ds, ps) -> int(round(n*ds[1]))),
   ((ds, ps) -> int(max(1.0, ceil(ps[1] / ds[2])))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4])
  ],
  [( 0.1, float(n)),
   ( 6.0, float(n)),
   (-2.0, 0.0),
   (0.01, 0.50)
  ]
)

outfile = "cmsa_es_exp8.csv"
write_csv_header_to_file(pe8, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 4 5 10 4 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe8, outfile; num_repeats = 1)
# Since many runs end up in the 2.0 magnitude class we optimize on fitness rather than magnitude => higher resolution.
# First we minimize the magnitude class with ei => 2*4+1
for i in 1:20
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 4 1 9 4 $(outfile) new_runs.json sa ei`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe8, outfile; num_repeats = 1)
end

# It seems hard to make progress on this one and it seems that very low values for lambda
# and mu are needed. So lets expand the search ranges and add back the old tau_c and tau
# covariates and lets see what we find.

pe9 = ParameterExperiment(
  ["lambda", "mu", "sigma", "decompose_covar_prob", "tau", "tau_c"
  ],
  [((ds, ps) -> int(10^ds[1])),
   ((ds, ps) -> int(max(1.0, ds[2]*ps[1]))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4]),
   ((ds, ps) -> (2*n)^ds[5]),
   ((ds, ps) -> 1 + n * (n+1) / (ds[6]*ps[2]))
  ],
  [( log10(4.0 + log(n)), log10(n*n)),
   ( 0.0, 1/4),
   (-3.0, 0.0),
   (0.01, 0.50),
   (-2.0, 0.0),
   (1.0,  6.0)
  ]
)
outfile = "cmsa_es_exp9.csv"
write_csv_header_to_file(pe9, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 6 7 14 6 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe9, outfile; num_repeats = 1)
for i in 1:23
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 6 1 14 6 $(outfile) new_runs.json sa ei`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe9, outfile; num_repeats = 1)
end
# Then we maximize the success rate and select the quickest among the top 10.
for i in 1:3
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 6 1 15 6 $(outfile) new_runs.json sa min`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe9, outfile; num_repeats = 1)
end

pe10 = ParameterExperiment(
  ["lambda", "mu", "sigma", "decompose_covar_prob",
   "max_successes_before_increasing", "max_failures_before_decreasing", "max_rounds_without_improvement"
  ],
  [((ds, ps) -> int(10^ds[1])),
   ((ds, ps) -> int(max(1.0, ds[2]*ps[1]))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4]),
   ((ds, ps) -> int(ds[5])),
   ((ds, ps) -> int(ds[6])),
   ((ds, ps) -> int(ds[7]))
  ],
  [( log10(4.0 + log(n)), log10(4*n)),
   ( 0.0, 1/7),
   (-2.0, 0.0),
   (0.01, 0.40),
   (2.0, 10.0),
   (2.0, 10.0),
   (20.0, 1000.0)
  ]
)
outfile = "cmsa_es_exp10.csv"
write_csv_header_to_file(pe10, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 8 16 7 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe10, outfile; num_repeats = 1)
for i in 1:13
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 16 7 $(outfile) new_runs.json sa ei`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe10, outfile; num_repeats = 1)
end
# Then we maximize the success rate and select the quickest among the top 10.
for i in 1:5
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 16 7 $(outfile) new_runs.json sa minquick`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe10, outfile; num_repeats = 1)
end

# Ok, we have solved Rosenbrock from 8-d to 128-d. Lets address something harder. Griewank?
n = 2
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Griewank"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe11 = ParameterExperiment(
  ["lambda", "mu", "sigma", "decompose_covar_prob",
   "max_successes_before_increasing", "max_failures_before_decreasing"
  ],
  [((ds, ps) -> int(10^ds[1])),
   ((ds, ps) -> int(max(1.0, ds[2]*ps[1]))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4]),
   ((ds, ps) -> int(ds[5])),
   ((ds, ps) -> int(ds[6]))
  ],
  [( log10(4.0 + log(n)), log10(8*n)),
   ( 1/30, 1/6),
   (-2.0, 0.0),
   (0.01, 0.99),
   (2.0, 1000.0),
   (2.0, 1000.0)
  ]
)
outfile = "cmsa_es_exp11.csv"
write_csv_header_to_file(pe11, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 6 7 14 6 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe11, outfile; num_repeats = 1)
for i in 1:13
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 6 1 13 6 $(outfile) new_runs.json sa ei`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe11, outfile; num_repeats = 1)
end

# Way too easy with 2-d Griewank though. Lets do 4-D and add all 7 params back in since the default value
# for max rounds without improvement sometimes terminates runs too early if max_failures_before_decreasing
# is also low.
n = 4
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Griewank"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe12 = ParameterExperiment(
  ["lambda", "mu", "sigma", "decompose_covar_prob",
   "max_successes_before_increasing", "max_failures_before_decreasing", "max_rounds_without_improvement"
  ],
  [((ds, ps) -> int(10^ds[1])),
   ((ds, ps) -> int(max(1.0, ds[2]*ps[1]))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4]),
   ((ds, ps) -> int(ds[5])),
   ((ds, ps) -> int(ds[6])),
   ((ds, ps) -> int(10^ds[7]))
  ],
  [( log10(4.0 + log(n)), log10(4*n*n)),
   ( 1/50, 1/2),
   (-2.0, 0.0),
   (0.01, 0.99),
   (-1000.0, 1000.0), # Inactive if less than 1
   (-1000.0, 1000.0), # Inactive if less than 1
   (2.0, 5.0)
  ]
)
outfile = "cmsa_es_exp12.csv"
write_csv_header_to_file(pe12, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 8 15 7 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe12, outfile; num_repeats = 1)
for i in 1:5
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 15 7 $(outfile) new_runs.json sa alc`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe12, outfile; num_repeats = 1)
end
for i in 1:10
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 16 7 $(outfile) new_runs.json sa min`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe12, outfile; num_repeats = 1)
end
for i in 1:5
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 15 7 $(outfile) new_runs.json sa alc`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe12, outfile; num_repeats = 1)
end

# The CMSA-ES performs rather badly here. Only one success in many tries, probably just by luck.
# Will need to try with some other CMA versions to see if they can manage with Griewank.

# Lets try a larger Griewank just to be sure it is generally bad for this.
n = 16
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Griewank"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe13 = ParameterExperiment(
  ["lambda", "mu", "sigma", "decompose_covar_prob",
   "max_successes_before_increasing", "max_failures_before_decreasing", "max_rounds_without_improvement"
  ],
  [((ds, ps) -> int(10^ds[1])),
   ((ds, ps) -> int(max(1.0, ds[2]*ps[1]))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4]),
   ((ds, ps) -> int(ds[5])),
   ((ds, ps) -> int(ds[6])),
   ((ds, ps) -> int(10^ds[7]))
  ],
  [( log10(4.0 + log(n)), log10(4*n*n)),
   ( 1/50, 1/2),
   (-3.0, 0.0),
   (0.01, 0.99),
   (-1000.0, 1000.0), # Inactive if less than 1
   (-1000.0, 1000.0), # Inactive if less than 1
   (2.0, 5.0)
  ]
)
outfile = "cmsa_es_exp13.csv"
write_csv_header_to_file(pe13, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 8 15 7 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe13, outfile; num_repeats = 10)
for i in 1:12
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 15 7 $(outfile) new_runs.json sa alc`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe13, outfile; num_repeats = 10)
end
for i in 1:10
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 19 7 $(outfile) new_runs.json sa alc`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe13, outfile; num_repeats = 10)
end
for i in 1:10
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 -19 7 $(outfile) new_runs.json sa minquick`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe13, outfile; num_repeats = 10)
end

# Strangely enough it works just fine for 16-D. Very strange. What is the problem with low-dim cases now; have
# I introduced some bias when focusing on the high-dim cases?

# Let's try 8D

n = 8
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Griewank"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe = pe14 = ParameterExperiment(
  ["lambda", "mu", "sigma", "decompose_covar_prob",
   "max_successes_before_increasing", "max_failures_before_decreasing", "max_rounds_without_improvement"
  ],
  [((ds, ps) -> int(10^ds[1])),
   ((ds, ps) -> int(max(1.0, ds[2]*ps[1]))),
   ((ds, ps) -> diameter * 10^ds[3]),
   ((ds, ps) -> ds[4]),
   ((ds, ps) -> int(ds[5])),
   ((ds, ps) -> int(ds[6])),
   ((ds, ps) -> int(10^ds[7]))
  ],
  [( log10(4.0 + log(n)), log10(4*n*n)),
   ( 1/50, 1/2),
   (-3.0, 0.0),
   (0.01, 0.99),
   (-1000.0, 1000.0), # Inactive if less than 1
   (-1000.0, 1000.0), # Inactive if less than 1
   (2.0, 5.0)
  ]
)
outfile = "cmsa_es_exp14.csv"
write_csv_header_to_file(pe14, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 8 15 7 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe, outfile; num_repeats = 10)
for i in 1:12
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 -19 7 $(outfile) new_runs.json not alc`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe, outfile; num_repeats = 10)
end
for i in 1:5
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 -19 7 $(outfile) new_runs.json sa minquick`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe, outfile; num_repeats = 5)
end

# It worked fine. Strange. lets try 4-d again.
n = 4
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Griewank"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

outfile = "cmsa_es_exp15.csv"
write_csv_header_to_file(pe, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 8 15 7 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe, outfile; num_repeats = 10)
for i in 1:12
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 -19 7 $(outfile) new_runs.json not alc`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe, outfile; num_repeats = 10)
end
for i in 1:5
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 -19 7 $(outfile) new_runs.json sa minquick`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe, outfile; num_repeats = 5)
end

# Now it works at least some of the time. So key was to do more than 1 run per setting... Lets try 2-D
n = 2
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Griewank"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

outfile = "cmsa_es_exp16.csv"
write_csv_header_to_file(pe, outfile)
run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 8 15 7 $(outfile) new_runs.json`)
run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe, outfile; num_repeats = 10)
for i in 1:12
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 -19 7 $(outfile) new_runs.json not alc`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe, outfile; num_repeats = 10)
end
for i in 1:5
  run(`/usr/bin/Rscript ../../R/parameter_experiment.R 7 1 -19 7 $(outfile) new_runs.json sa minquick`)
  run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(cmsa_es, p, pe, outfile; num_repeats = 5)
end

# Looks better. Ok, lets try Sphere 2-D with the new experiment function that
# does it all. Just the simplest test to make sure everything works.

n = 2
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Sphere"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

pe = pe15 = ParameterExperiment(
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
explore_parameters(pe15, cmsa_es, p; experiment_prefix = "cmsa_es", num_repeats = 10)
