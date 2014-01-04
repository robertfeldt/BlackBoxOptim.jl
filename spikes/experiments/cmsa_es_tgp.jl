using BlackBoxOptim
using JSON

require("../covar_matrix_self_adaptation_es.jl")
require("../experiment_framework.jl")

type ParameterExperiment
  parameters::Array{ASCIIString, 1}           # Parameter names
  mappers::Array{Function, 1}                 # One function per param mapping a design value to a param value
  design_ranges::Array{(Float64, Float64), 1}   # Design range per parameter
end

n = 8
p = as_fixed_dim_problem(BlackBoxOptim.example_problems["Rosenbrock"], n)
p = BlackBoxOptim.shifted(p)
diameter = minimum(diameters(search_space(p)))

cmsa_es_exp1 = ParameterExperiment(
  ["lambda", "mu", "covarMatrixSampler", "utilitiesFunc", 
   "covar_learning_rate", "sigma"],
  [((ds, ps) -> 4 * [2, n, n*n][int(ceil(ds[1]))]),
   ((ds, ps) -> int(ps[1] / [2, 4, 8][int(ceil(ds[2]))])),
   ((ds, ps) -> [EigenCovarSampler, CholeskyCovarSampler][int(ceil(ds[3]))]),
   ((ds, ps) -> [log_utilities, linear_utilities][int(ceil(ds[4]))]),
   ((ds, ps) -> ds[5]),
   ((ds, ps) -> diameter * 10^ds[6])
  ],
  [( 0.0, 3.0),
   ( 0.0, 3.0),
   ( 0.0, 2.0),
   ( 0.0, 2.0),
   ( 0.6, 1.0),
   (-2.0, 0.0)]
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

function design_ranges_to_R_array(pe::ParameterExperiment)
  ary = Float64[]
  for(i in 1:numparams(pe))
    if !is_fixed_param(pe, i)
      push!(ary, float(pe.design_ranges[i][1]))
      push!(ary, float(pe.design_ranges[i][2]))
    end
  end
  array_to_R_string(ary)
end

function array_to_R_string(a)
  join(["c(", join(a, ","), ")"])
end

function write_csv_header_to_file(pe::ParameterExperiment, filepath)
  header = join([
    map((i) -> "d$(i)", 1:numparams(pe)), 
    pe.parameters,
    ["Fitness", "NumFuncEvals", "Time", "Reason"]],
    ",")
  fh = open(filepath, "w")
  println(fh, header)
  close(fh)
end

# Now we want to invoke an R script like so:
# Rscript top_of_bbo_path/R/parameter_experiment.R 6 "c(0.0,3.0,0.0,3.0,0.0,2.0,0.0,2.0,0.6,1.0,-2.0,0.0)" 7 13 in.csv out.json
# with parameters indicating: 6 = num params, ranges, num runs, response col, in.csv, out.json
# where in.csv has the values used so far and out.json has the 7*6 new design values.

function read_matrix_from_file(filename, nrows, ncols)
  fh = open(filename, "r")
  fromR = JSON.parse(readall(fh))
  close(fh)
  reshape(fromR, nrows, ncols)
end

function create_design_matrix_from_file_and_fixed_params(pe::ParameterExperiment, filename, num_rows)
  design = read_matrix_from_file(filename, num_rows, num_varying_params(pe))

  ds = zeros(num_rows, numparams(pe))
  next_design_col = 1

  for(i in 1:numparams(pe))
    if is_fixed_param(pe, i)
      ds[:, i] = value_for_fixed_param(pe, i) * ones(num_rows)
    else
      ds[:, i] = design[:,next_design_col]
      next_design_col += 1
    end
  end

  ds
end

function run_based_on_design_matrix_while_saving_to_csvfile(pe::ParameterExperiment, design, outfilepath)

  for(i in 1:size(design,1))
    ds = design[i,:]
    ps = Any[]
    param_dict = Dict{Any,Any}()

    # Calc the parameter values from the desing values.
    for(j in 1:numparams(pe))
      param_dict[pe.parameters[j]] = pv = pe.mappers[j](ds, ps)
      push!(ps, pv)
    end

    # Add other arguments.
    param_dict[:max_evals_per_dim] = 1e7

    # Now run it while timing.
    tic()
    xb, fb, nf, r, a = cmsa_es(p; collect(param_dict)...)
    t = toq()

    # Save info to the csv file
    fh = open(outfilepath, "a+")
    print(fh, join(map( (dv) -> "$(dv)", ds), ","))
    print(fh, ",")
    print(fh, join(map( (pv) -> "$(pv)", ps), ","))
    println(fh, ",", fb, ",", nf, ",", t, ",", r)
  end

end
