"""
Create `Evaluator` instance for a given `problem`.
"""
function make_evaluator(problem::OptimizationProblem, archive=nothing, params::Parameters=ParamsDict())
    workers = get(params, :Workers, Vector{Int}())
    if archive===nothing
        # make the default archive
        archiveCapacity = get(params, :ArchiveCapacity, 10)
        archive = TopListArchive(fitness_scheme(problem), numdims(problem), archiveCapacity)
    end
    if length(workers) > 0
        return ParallelEvaluator(problem, archive, pids=workers)
    else
        return ProblemEvaluator(problem, archive)
    end
end

"""
Optimization Run Controller.
Manages problem optimization using the specified method.

See `OptController`.
"""
mutable struct OptRunController{O<:Optimizer, E<:Evaluator}
    optimizer::O   # optimization algorithm
    evaluator::E   # problem evaluator

    trace_mode::Symbol # controller state trace mode (:verbose, :compact, :silent)
                       # :silent makes tr() generate no output)
    save_trace::Bool # FIXME if traces should be saved to a file
    trace_interval::Float64 # periodicity of calling trace_progress()

    callback_function::Function
    callback_interval::Float64

    # termination conditions
    max_steps::Int      # maximal number of steps
    max_fevals::Int # maximal number of fitness evaluations
    max_steps_without_fevals::Int # stop optimization if no func evals in this many steps (indicates a converged/degenerate search)
    max_steps_without_progress::Int # stop optimization if no improvement in this many steps (indicates a converged/degenerate search)
    max_time::Float64   # maximal time, 0 to ignore

    min_delta_fitness_tol::Float64 # minimal difference between current best fitness and second-best one
    fitness_tol::Float64  # minimal difference between current best fitness and the known optmimum

    num_steps::Int       # optimization steps done
    num_better::Int      # number of steps that improved best fitness

    num_better_since_last_report::Int # ditto, since last call to trace_progress()
    num_steps_since_last_report::Int # number of steps since last call to trace_progress()

    last_num_fevals::Int # the number of function evals on the previous step
    num_steps_without_fevals::Int # the number of steps without the function evals

    start_time::Float64 # time optimization started, 0 if not running yet
    stop_time::Float64  # time optimization stopped, 0 if still running
    last_report_time::Float64 # last time trace_progress() was called
    last_callback_time::Float64 # last time callback function was called

    stop_reason::String # the reason for algorithm termination, empty if it's not terminated
end

"""
Create `OptRunController` for a given problem using specified `optimizer`.

#Arguments
    * `optimizer` initialized optimization method
    * `evaluator` the evaluator of the problem fitness
    * `params` controller settings, see `DefaultParameters` for the default values:
        * `:MaxTime` max time in seconds (takes precedence over the other budget-related params if specified), 0.0 disables the check
        * `:MaxFuncEvals` max fitness evals (takes precedence over max iterations, but not max time), 0 disables the check
        * `:MaxSteps` max iterations gives the least control since different optimizers have different "size" of their "iterations"
        * `:MaxStepsWithoutProgress` max iterations without fitness improvement
        * `:MinDeltaFitnessTolerance` minimum delta fitness (difference between the two consecutive best fitness improvements) we can accept before terminating
        * `:FitnessTolerance` stop the optimization when the best fitness found is within this distance of the actual optimum (if known)
        * `:MaxNumStepsWithoutFuncEvals` stop optimization if no new fitness evals in this many steps (indicates a converged/degenerate search)
        * `:NumRepetitions` number of repetitions to run for each optimizer for each problem
        * `:TraceMode` how the optimizer state is traced to the STDOUT during the optimization (one of `:silent`, `:verbose`)
        * `:TraceInterval` the trace interval (in seconds)
        * `:SaveTrace` whether to save it to a file (defaults to `false`)
        * `:SaveFitnessTraceToCsv` whether the history of fitness changes during optimization should be save to a csv file
        * `:SaveParameters` save method/controller parameters to a JSON file
"""
function OptRunController(optimizer::O, evaluator::E, params) where {O<:Optimizer, E<:Evaluator}
    OptRunController{O,E}(optimizer, evaluator,
        [params[key] for key in Symbol[:TraceMode, :SaveTrace, :TraceInterval,
                      :CallbackFunction, :CallbackInterval,
                      :MaxSteps, :MaxFuncEvals, :MaxNumStepsWithoutFuncEvals, :MaxStepsWithoutProgress, :MaxTime,
                      :MinDeltaFitnessTolerance, :FitnessTolerance]]...,
        0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, -1.0, "")
end

# stepping optimizer has it's own evaluator, get a reference
OptRunController(optimizer::SteppingOptimizer, problem::OptimizationProblem, params) =
    OptRunController(optimizer, evaluator(optimizer), params)
# all other optimizers are using make_evaluator() method to create evaluator by default
OptRunController(optimizer::Optimizer, problem::OptimizationProblem, params) =
    OptRunController(optimizer, make_evaluator(problem, nothing, params), params)
# We reuse evaluator and optimizer from a previous run if set up from such a run.
OptRunController(orc::OptRunController, params) =
    OptRunController(optimizer(orc), evaluator(orc), params)

# logging/tracing
function trace(ctrl::OptRunController, msg::AbstractString, obj = nothing)
    if ctrl.trace_mode != :silent
        print(msg)
        if obj !== nothing
            if isa(obj, AbstractString)
                print(obj)
            else
                show(IOContext(stdout, :compact => true), obj)
            end
        end
    end
    if ctrl.save_trace
        # No saving for now
    end
end

optimizer(ctrl::OptRunController) = ctrl.optimizer
evaluator(ctrl::OptRunController) = ctrl.evaluator
problem(ctrl::OptRunController) = problem(evaluator(ctrl))
isstarted(ctrl::OptRunController) = ctrl.start_time > 0
isrunning(ctrl::OptRunController) = isstarted(ctrl) && ctrl.stop_time == 0
isstopped(ctrl::OptRunController) = ctrl.stop_time > 0

num_steps(ctrl::OptRunController) = ctrl.num_steps
num_func_evals(ctrl::OptRunController) = num_evals(ctrl.evaluator)
stop_reason(ctrl::OptRunController) = ctrl.stop_reason
start_time(ctrl::OptRunController) = ctrl.start_time

best_candidate(ctrl::OptRunController) = best_candidate(ctrl.evaluator.archive)
best_fitness(ctrl::OptRunController) =  best_fitness(ctrl.evaluator.archive)

"""
    show_fitness(io, fit, [problem::OptimizationProblem])

Output fitness to the given IO stream.
`show_fitness()` method could be overridden for a specific problem, e.g.
to print the names of each objective.
"""
show_fitness(io::IO, fit::Number) = @printf(io, "%.9f", fit)

function show_fitness(io::IO, fit::NTuple{N}) where N
    print(io, "(")
    for i in 1:N
        if i > 1
            print(io, ", ")
        end
        @printf(io, "%.5f", fit[i])
    end
    print(io, ")")
end

function show_fitness(io::IO, fit::IndexedTupleFitness{N}) where N
    show_fitness(io, fit.orig)
    @printf(io, " agg=%.5f", fit.agg)
end

# problem-specific method defaults to problem-agnostic method
# note that the fitness type of the problem could be different from FA,
# because show_fitness() typically is called for the archived fitness
show_fitness(io::IO, fit::FA, problem::OptimizationProblem) where {FA} = show_fitness(io, fit)

"""
    format_fitness(fit, [problem::OptimizationProblem])

Format fitness into a string.
Calls `show_fitness()` under the hood.
"""
function format_fitness(fit::Any, problem::OptimizationProblem)
    buf = IOBuffer(false, true)
    show_fitness(buf, fit, problem)
    return bytestring(buf.data)
end

function format_fitness(fit::Any)
    buf = IOBuffer(false, true)
    show_fitness(buf, fit)
    return bytestring(buf.data)
end

function elapsed_time(ctrl::OptRunController)
    isrunning(ctrl) ? time() - ctrl.start_time : ctrl.stop_time - ctrl.start_time
end

function check_stop_condition(ctrl::OptRunController)
    if ctrl.max_time > 0 && elapsed_time(ctrl) > ctrl.max_time
        return "Max time ($(ctrl.max_time) s) reached"
    end

    if ctrl.max_fevals > 0 && num_func_evals(ctrl) > ctrl.max_fevals
        return "Max number of function evaluations ($(ctrl.max_fevals)) reached"
    end

    if ctrl.max_steps_without_fevals > 0 && ctrl.num_steps_without_fevals > ctrl.max_steps_without_fevals
        return "Too many steps ($(ctrl.num_steps_without_fevals)) without any function evaluations (probably search has converged)"
    end

    if ctrl.max_steps > 0 && num_steps(ctrl) > ctrl.max_steps
        return "Max number of steps ($(ctrl.max_steps)) reached"
    end

    return check_stop_condition(ctrl.evaluator, ctrl)

    return "" # empty string, no termination
end

@inline function callback(ctrl::OptRunController)
    if ctrl.callback_interval > 0.0
        ctrl.callback_function(ctrl)
        ctrl.last_callback_time = time()
    end
end

function trace_progress(ctrl::OptRunController)
    # update the counters
    ctrl.last_report_time = time()
    total_improvement_rate = ctrl.num_better/num_steps(ctrl)
    recent_improvement_rate = ctrl.num_better_since_last_report/ctrl.num_steps_since_last_report
    ctrl.num_better_since_last_report = 0
    ctrl.num_steps_since_last_report = 0

    if ctrl.trace_mode == :silent
        return
    end

    # Always print step number, num fevals and elapsed time
    trace(ctrl, @sprintf("%.2f secs, %d evals, %d steps",
            elapsed_time(ctrl), num_func_evals(ctrl), num_steps(ctrl)))

    # Only print if this optimizer reports on number of better. They return 0
    # if they do not.
    if total_improvement_rate > 0.0
        trace(ctrl, @sprintf(", improv/step: %.3f (last = %.4f)",
                total_improvement_rate, recent_improvement_rate))
    end

    # Always print fitness if num_evals > 0
    if num_func_evals(ctrl) > 0
        trace(ctrl, ", fitness=")
        show_fitness(stdout, best_fitness(ctrl), problem(ctrl))
    end

    trace(ctrl, "\n")

    trace_state(stdout, ctrl.optimizer, ctrl.trace_mode)
end

function step!(ctrl::OptRunController{<:AskTellOptimizer})
    # The ask()/tell() interface is more general since you can mix and max
    # elements from several optimizers using it. However, in this top-level
    # execution function we do not make use of this flexibility...
    candidates = ask(ctrl.optimizer)
    rank_by_fitness!(ctrl.evaluator, candidates)
    return tell!(ctrl.optimizer, candidates)
end

# step for SteppingOptimizers
function step!(ctrl::OptRunController{<:SteppingOptimizer})
    step!(ctrl.optimizer)
    return 0
end

setup_optimizer!(ctrl::OptRunController{<:SteppingOptimizer}) =
    setup!(ctrl.optimizer)
setup_optimizer!(ctrl::OptRunController{<:AskTellOptimizer}) =
    setup!(ctrl.optimizer, ctrl.evaluator)

shutdown_optimizer!(ctrl::OptRunController{<:SteppingOptimizer}) =
    shutdown!(ctrl.optimizer)

function shutdown_optimizer!(ctrl::OptRunController{<:AskTellOptimizer})
    shutdown!(ctrl.optimizer)
    shutdown!(ctrl.evaluator)
end

"""
    run!(ctrl::OptRunController)

Run optimization until one of the stopping conditions are satisfied.
"""
function run!(ctrl::OptRunController)
    trace(ctrl, "Starting optimization with optimizer $(name(ctrl.optimizer))\n")
    try
        setup_optimizer!(ctrl)

        ctrl.start_time = time()
        ctrl.num_steps = 0
        while isempty(ctrl.stop_reason)
            # Report on progress every now and then...
            if (time() - ctrl.last_report_time) > ctrl.trace_interval
                trace_progress(ctrl)
            end

            # Take the step and then update the counters
            nstep_better = step!(ctrl)
            ctrl.num_better += nstep_better
            ctrl.num_better_since_last_report += nstep_better
            ctrl.num_steps += 1
            ctrl.num_steps_since_last_report += 1
            if num_evals(ctrl.evaluator) == ctrl.last_num_fevals
                ctrl.num_steps_without_fevals += 1
            else
                ctrl.num_steps_without_fevals = 0
            end
            ctrl.last_num_fevals = num_func_evals(ctrl)

            # Callback every now and then (if a callback interval has been set)...
            if ctrl.callback_interval > 0.0 && 
                (ctrl.last_callback_time <= 0.0 ||
                    (time() - ctrl.last_callback_time) > ctrl.callback_interval)
                callback(ctrl)
            end            
            
            # Check if we have a reason to stop on next loop
            ctrl.stop_reason = check_stop_condition(ctrl)
        end
        ctrl.stop_time = time()
    finally
        shutdown_optimizer!(ctrl)
    end
    trace(ctrl, "\nOptimization stopped after $(ctrl.num_steps) steps and $(elapsed_time(ctrl)) seconds\n")

    return ctrl.stop_reason
end

function show_report(ctrl::OptRunController, population_stats=false)
    final_elapsed_time = elapsed_time(ctrl)
    trace(ctrl, "Termination reason: $(ctrl.stop_reason)\n")
    trace(ctrl, "Steps per second = $(num_steps(ctrl)/final_elapsed_time)\n")
    trace(ctrl, "Function evals per second = $(num_func_evals(ctrl)/final_elapsed_time)\n")
    trace(ctrl, "Improvements/step = $(ctrl.num_better/ctrl.max_steps)\n")
    trace(ctrl, "Total function evaluations = $(num_func_evals(ctrl))\n")

    if population_stats && isa(ctrl.optimizer, PopulationOptimizer)
        trace(ctrl, "\nMean value (in population) per position:",  mean(population(ctrl.optimizer),1))
        trace(ctrl, "\n\nStd dev (in population) per position:", std(population(ctrl.optimizer),1))
    end

    trace(ctrl, "\n\nBest candidate found: ", best_candidate(ctrl))
    trace(ctrl, "\n\nFitness: ")
    show_fitness(stdout, best_fitness(ctrl), problem(ctrl))
    trace(ctrl, "\n\n")
end

function write_result(ctrl::OptRunController, filename = "")
    if isempty(filename)
        timestamp = strftime("%y%m%d_%H%M%S", floor(Int, ctrl.start_time))
        filename = "$(timestamp)_$(problem_summary(ctrl.evaluator))_$(name(ctrl.optimizer)).csv"
        filename = replace(replace(filename, r"\s+", "_"), r"/", "_")
    end
    save_fitness_history_to_csv_file(ctrl.evaluator.archive, filename;
            header_prefix = "Problem,Dimension,Optimizer",
            line_prefix = "$(name(problem(ctrl.evaluator))),$(numdims(ctrl.evaluator)),$(name(ctrl.optimizer))",
            bestfitness = opt_value(problem(ctrl.evaluator)))
end

"""
Optimization Controller.

Applies specific optimization method to a given problem.
Supports restarts and modifying parameter of the method between runs.
`runcontrollers` field maintains the list of `OptRunController` instances applied so far.

See `OptRunController`.
"""
mutable struct OptController{O<:Optimizer, P<:OptimizationProblem}
    optimizer::O   # optimization algorithm
    problem::P     # opt problem
    parameters::ParamsDictChain
    runcontrollers::Vector{OptRunController{O}}
end

"""
Create `OptController` for a given `optimizer` and a `problem`.

`params` are any of `OptRunController` parameters plus
    * `:RngSeed` and `:RandomizeRngSeed` params for controlling the random seed
    * `:RecoverResults` if intermediate results are returned upon `InterruptException()` (on by default)
"""
OptController(
    optimizer::O, problem::P,
    params::ParamsDictChain) where {O<:Optimizer, P<:OptimizationProblem} =
    OptController{O, P}(optimizer, problem, params, OptRunController{O}[])

problem(oc::OptController) = oc.problem
numruns(oc::OptController) = length(oc.runcontrollers)
lastrun(oc::OptController) = oc.runcontrollers[end]

"""
    update_parameters!(oc::OptController, parameters::Associative)

Update the `OptController` parameters.
"""
function update_parameters!(oc::OptController, parameters::Parameters = EMPTY_DICT)
    parameters = convert(ParamsDict, parameters)

    # Most parameters cannot be changed since the problem and optimizer has already
    # been setup.
    valid_params = Set([:MaxTime, :MaxSteps, :MaxFuncEvals, :TraceMode])
    for k in keys(parameters)
        if k ∉ valid_params
            throw(ArgumentError("It is currently not supported to change parameters that can affect the original opt problem or optimizer (here: $(k))"))
        end
    end

    # Add new params in front if any are specified and they are valid.
    if length(parameters) > 0
        oc.parameters = chain(oc.parameters, parameters)
        check_valid!(oc.parameters) # We must recheck that new param settings are valid
    end
end

function init_rng!(parameters::Parameters)
    if parameters[:RandomizeRngSeed]
        parameters[:RngSeed] = rand(1:1_000_000)
        Random.seed!(parameters[:RngSeed])
    end
end

"""
    run!(oc::OptController)

Start a new optimization run, possibly with new parameters and report on results.
"""
function run!(oc::OptController)
    # If this is not the first run we create the controller based on the previous one
    # to reuse the archive etc, otherwise create it based on the problem.
    ctrl = if numruns(oc) > 0
        OptRunController(lastrun(oc), oc.parameters)
    else
        # No previous run controller so create a new one from optimizer and problem.
        OptRunController(oc.optimizer, oc.problem, oc.parameters)
    end
    push!(oc.runcontrollers, ctrl)

    # If this is the first run we might have to init the RNG.
    if numruns(oc) == 1
        init_rng!(oc.parameters)
    end

    try
        run!(ctrl)
        if ctrl.trace_mode != :silent
            show_report(ctrl)
        end

        if oc.parameters[:SaveFitnessTraceToCsv]
            write_results(ctrl)
        end

        return OptimizationResults(ctrl, oc)
    catch ex
        # If it was a ctrl-c interrupt we try to make a result and return it...
        if get(oc.parameters, :RecoverResults, true) && isa(ex, InterruptException)
            @warn("Optimization interrupted, recovering intermediate results...")
            ctrl.stop_reason = @sprintf("%s", ex)
            return OptimizationResults(ctrl, oc)
        else
            rethrow(ex)
        end
    end
end
