using BlackBoxOptim
using JSON

require("../covar_matrix_self_adaptation_es.jl")
require("../experiment_framework.jl")

type ParameterExperiment
  parameters::Array{ASCIIString, 1}             # Parameter names
  mappers::Array{Function, 1}                   # One function per param mapping a design value to a param value
  design_ranges::Array{(Float64, Float64), 1}   # Design range per parameter
end

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

# Return true iff a param with a given index is fixed, i.e. if its range
# has same min as max.
function is_fixed_param(pe::ParameterExperiment, index::Int64)
  pe.design_ranges[index][1] == pe.design_ranges[index][2]
end

numparams(pe::ParameterExperiment) = length(pe.parameters)
num_varying_params(pe::ParameterExperiment) = length(index_to_varying_params(pe))
num_fixed_params(pe::ParameterExperiment) = numparams(pe) - num_varying_params(pe)
index_to_varying_params(pe::ParameterExperiment) = filter( (i) -> !is_fixed_param(pe, i), 1:numparams(pe) )
value_for_fixed_param(pe::ParameterExperiment, i) = pe.design_ranges[i][1]

#function design_ranges_to_R_array(pe::ParameterExperiment)
#  ary = Float64[]
#  for(i in 1:numparams(pe))
#    if !is_fixed_param(pe, i)
#      push!(ary, float(pe.design_ranges[i][1]))
#      push!(ary, float(pe.design_ranges[i][2]))
#    end
#  end
#  array_to_R_string(ary)
#end
#
#function array_to_R_string(a)
#  join(["c(", join(a, ","), ")"])
#end

function write_csv_header_to_file(pe::ParameterExperiment, filepath)
  header = join([
    map((i) -> "d$(i)", 1:numparams(pe)), 
    pe.parameters,
    ["MedianFitness", "MedianNumFuncEvals", "SuccessRate"]],
    ",")
  fh = open(filepath, "w")
  println(fh, header)
  close(fh)
end

# Now we want to invoke an R script like so:
# Rscript top_of_bbo_path/R/parameter_experiment.R 6 7 13 in.csv out.json
# with parameters indicating: 6 = num params, ranges, num runs, response col, in.csv, out.json
# where in.csv has the values used so far and out.json has the 7*6 new design values.

function read_matrix_from_file(filename, nrows, ncols)
  fh = open(filename, "r")
  fromR = JSON.parse(readall(fh))
  close(fh)
  reshape(fromR, nrows, ncols)
end

function map01_range_to_design_range(value01, design_min, design_max)
  design_min + value01 * (design_max - design_min)
end

function create_design_matrix_from_file_and_fixed_params(pe::ParameterExperiment, filename, num_rows)
  design01 = read_matrix_from_file(filename, num_rows, num_varying_params(pe))

  ds = zeros(num_rows, numparams(pe))
  next_design_col = 1

  for(i in 1:numparams(pe))
    if is_fixed_param(pe, i)
      ds[:, i] = value_for_fixed_param(pe, i) * ones(num_rows)
    else
      dmin, dmax = pe.design_ranges[i]
      ds[:, i] = map01_range_to_design_range(design01[:,next_design_col], dmin, dmax)
      next_design_col += 1
    end
  end

  return design01, ds
end

function run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(pe::ParameterExperiment, outfilepath, num_rows; 
  designfile = "new_runs.json", num_repeats = 10, 
  fixed_params = {:max_evals_per_dim => 1e7, 
    :covarMatrixSampler => CholeskyCovarSampler,
    :utilitiesFunc => log_utilities})

  if !isfile(outfilepath)
    write_csv_header_to_file(pe, outfilepath)
  end

  design01, design = create_design_matrix_from_file_and_fixed_params(pe, designfile, num_rows)

  for(i in 1:size(design,1))
    ds = design[i,:]
    ps = Any[]
    param_dict = Dict{Any,Any}()

    # Calc the parameter values from the desing values.
    for(j in 1:numparams(pe))
      param_dict[symbol(pe.parameters[j])] = pv = pe.mappers[j](ds, ps)
      push!(ps, pv)
    end

    # Add other arguments.
    param_dict = merge(param_dict, fixed_params)
    print("params: "); show(param_dict); println("")

    # Set up for saving results
    fbs = zeros(num_repeats)
    nfs = zeros(num_repeats)
    num_within_ftol = 0

    # Now run it while timing.
    for(r in 1:num_repeats)
      tic()
      xb, fbs[r], nfs[r], r, a = cmsa_es(p; collect(param_dict)...)
      et = toq()
      if r == "Within ftol"
        num_within_ftol += 1
      end
    end

    # Save info to the csv file
    fh = open(outfilepath, "a+")
    print(fh, join(map( (dv) -> "$(dv)", design01[i,:]), ","))
    print(fh, ",")
    print(fh, join(map( (pv) -> "$(pv)", ps), ","))
    println(fh, ",", median(fbs), ",", median(nfs), ",", num_within_ftol/num_repeats)
    close(fh)
  end

end

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
