using JSON

struct ParameterExperiment
    parameters::Vector{String}               # Parameter names
    mappers::Vector{Function}                     # One function per param mapping a design value to a param value
    design_ranges::Vector{Tuple{Float64, Float64}}# Design range per parameter
end

numparams(pe::ParameterExperiment) = length(pe.parameters)

function write_csv_header_to_file(pe::ParameterExperiment, filepath)
    header = join([
        map((i) -> "d$(i)", 1:numparams(pe)),
        pe.parameters,
        ["MedianFitness", "MedianFitnessMagnitude", "MedianNumFuncEvals", "MedianTime", "SuccessRate"]],
        ",")
    fh = open(filepath, "w")
    println(fh, header)
    close(fh)
end

# The result from running the R tgp script is saved as a list to a json file
# that we read into a dict.
function read_json_design_list_from_file(filename)
    fh = open(filename, "r")
    fromR = JSON.parse(readall(fh))
    close(fh)
    fromR
end

function design01_matrix_from_R_list(fromR, matrix_name = "best")
    reshape(fromR[matrix_name], fromR["num_rows"], fromR["num_cols"])
end

function map01_range_to_design_range(value01, design_min, design_max)
    design_min + value01 * (design_max - design_min)
end

function design_matrix_from_design01_matrix(pe::ParameterExperiment, design01)
    ds = zeros(Float64, size(design01))
    next_design_col = 1

    for i in 1:numparams(pe)
        #if is_fixed_param(pe, i)
        #  ds[:, i] = value_for_fixed_param(pe, i) * ones(num_rows)
        #else
        dmin, dmax = pe.design_ranges[i]
        ds[:, i] = map01_range_to_design_range(design01[:,next_design_col], dmin, dmax)
        next_design_col += 1
        #end
    end

    return ds
end

function parameter_values_from_design_matrix(pe::ParameterExperiment, design)
    pvs = Any[]

    for i in 1:size(design,1)
        ds = design[i,:]
        ps = Any[]
        param_dict = ParamsDict()

        # Calc the parameter values from the design values.
        for j in 1:numparams(pe)
            param_dict[symbol(pe.parameters[j])] = pv = pe.mappers[j](ds, ps)
            push!(ps, pv)
        end

        push!(pvs, param_dict)
    end

    return pvs
end

function summary_of_R_output(pe::ParameterExperiment, designfile = "new_runs.json")
    res = read_json_design_list_from_file(designfile)
    for d01name in ["best", "top5_max", "top5_mean", "top5_min"]
        d01 = design01_matrix_from_R_list(res, d01name)
        d = design_matrix_from_design01_matrix(pe, d01)
        pvs = parameter_values_from_design_matrix(pe, d)
        print(d01name, ": ", pvs)
    end
    perm = sortperm(res["sa_mean_1st_order_sens_indices"], rev=true)
    println("Sensitivity analysis, decreasing main effects: ", pe.parameters[perm])
    perm = sortperm(res["sa_mean_total_sens_indices"], rev=true)
    println("Sensitivity analysis, decreasing total effects: ", pe.parameters[perm])
end

function run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(
        runfunc,
        problem, pe::ParameterExperiment, outfilepath;
        designfile = "new_runs.json", d01name = "best", num_repeats = 5,
        fixed_params = {:max_evals_per_dim => 1e7,
                        :utilitiesFunc => log_utilities}
)
    if !isfile(outfilepath)
        write_csv_header_to_file(pe, outfilepath)
    end

    res = read_json_design_list_from_file(designfile)
    d01 = design01_matrix_from_R_list(res, d01name)
    d = design_matrix_from_design01_matrix(pe, d01)
    pvs = parameter_values_from_design_matrix(pe, d)

    for i in 1:length(pvs)
        param_dict = pvs[i]
        # Add other arguments.
        param_dict = chain(param_dict, fixed_params)
        print("params: "); show(param_dict); println("")

        # Set up for saving results
        fbs = zeros(num_repeats)
        nfs = zeros(num_repeats)
        ets = zeros(num_repeats)
        num_within_ftol = 0

        # Now run it while timing.
        for k in 1:num_repeats
            tic()
            xb, fbs[k], nfs[k], r, a = runfunc(problem; collect(param_dict)...)
            ets[k] = toq()
            if r == "Within ftol"
                num_within_ftol += 1
            end
        end

        # Save info to the csv file
        fh = open(outfilepath, "a+")
        print(fh, join(map( (dv) -> "$(dv)", d01[i,:]), ","))
        for i in 1:numparams(pe)
            pv = param_dict[symbol(pe.parameters[i])]
            print(fh, ",$(pv)")
        end
        mf = median(fbs)
        fmagnitude = BlackBoxOptim.magnitude_class(mf)[2]
        println(fh, ",", median(fbs), ",", fmagnitude, ",", median(nfs), ",", median(ets), ",", num_within_ftol/num_repeats)
        close(fh)
    end
end

function run_with_parameters_from_json_file(pe::ParameterExperiment, designfile, searchfunc, problem;
    outfile = "exp.csv", num_repeats = 1, path_to_R_dir = "../../R")

    if !isfile(outfile)
        write_csv_header_to_file(pe, outfile)
    end

    run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(searchfunc,
        problem, pe, outfile; num_repeats = num_repeats, designfile = designfile)
end

function explore_parameters(
        pe::ParameterExperiment, searchfunc, problem;
        experiment_prefix = "exp", num_repeats = 10,
        path_to_R_dir = "../../R"
)
    nps = numparams(pe)
    experiment_name = join([experiment_prefix, name(problem), numdims(problem),
            nps, "params"], "_")

    outfile = join([experiment_name, ".csv"])
    designfile = join(["new_runs_for_", experiment_name, ".csv"])

    write_csv_header_to_file(pe, outfile)

    index_fitness = 2*nps+1
    index_magnitude = 2*nps+2
    index_fevals = 2*nps+3
    index_time = 2*nps+4
    index_success_rate = 2*nps+5

    # First run nps+1 times with a random LHS sample of points. Before we have
    # that number of points we cannot run the GP regressions.
    run(`/usr/bin/Rscript $(path_to_R_dir)/parameter_experiment.R $(nps) $(nps+1) $(index_fitness) $(nps) $(outfile) $(designfile)`)
    run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(searchfunc,
            problem, pe, outfile; num_repeats = num_repeats, designfile = designfile)

    # Now run new points, each selected by ALC and optimizing fitness magnitude,
    # until we have gathered 30 points in total.
    for i in 1:(30-nps-1)
        run(`/usr/bin/Rscript $(path_to_R_dir)/parameter_experiment.R $(nps) 1 $(index_magnitude) $(nps) $(outfile) $(designfile) not alc`)
        run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(searchfunc,
                problem, pe, outfile; num_repeats = num_repeats, designfile = designfile)
    end

    # Now run 20 new points that minimize fitness and select by alc.
    for i in 1:20
        run(`/usr/bin/Rscript $(path_to_R_dir)/parameter_experiment.R $(nps) 1 $(index_fitness) $(nps) $(outfile) $(designfile) not alc`)
        run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(searchfunc,
                problem, pe, outfile; num_repeats = num_repeats, designfile = designfile)
    end

    # Now run 10 new points that minimize fitness magnitude and select greedily to minimize execution time.
    for i in 1:10
        run(`/usr/bin/Rscript $(path_to_R_dir)/parameter_experiment.R $(nps) 1 $(index_magnitude) $(nps) $(outfile) $(designfile) not min`)
        run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(searchfunc,
                problem, pe, outfile; num_repeats = num_repeats, designfile = designfile)
    end

    # Now run 10 new points that minimize fitness magnitude and select greedily to
    # minimize and be quick.
    for i in 1:10
        run(`/usr/bin/Rscript $(path_to_R_dir)/parameter_experiment.R $(nps) 1 $(index_magnitude) $(nps) $(outfile) $(designfile) not min`)
        run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(searchfunc,
                problem, pe, outfile; num_repeats = num_repeats, designfile = designfile)
    end

    # Now run 10 new points that maximize success rate and select greedily to
    # minimize and be quick. I missed the negative sign on the first runs so skip
    # those runs.
    for i in 1:10
        run(`/usr/bin/Rscript $(path_to_R_dir)/parameter_experiment.R $(nps) 1 -$(index_success_rate) $(nps) $(outfile) $(designfile) not minquick`)
        run_based_on_design_matrix_in_file_while_saving_results_to_csvfile(searchfunc,
                problem, pe, outfile; num_repeats = num_repeats, designfile = designfile)
    end
end
