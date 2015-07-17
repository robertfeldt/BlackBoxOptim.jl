# Optimization Controller
# Applies specific optimizer to a given optimization problem.
type OptController{O<:Optimizer, E<:Evaluator}
  optimizer::O   # optimization algorithm
  evaluator::E   # problem evaluator

  show_trace::Bool # if controller should trace its execution (false makes tr() generate no output)
  save_trace::Bool # FIXME if traces should be saved to a file
  trace_interval::Float64 # periodicity of calling trace_progress()

  # termination conditions
  max_steps::Int      # maximal number of steps
  max_fevals::Int # maximal number of fitness evaluations
  max_steps_without_fevals::Int # stop optimization if no func evals in this many steps (indicates a converged/degenerate search)
  max_time::Float64   # maximal time, 0 to ignore

  min_delta_fitness_tol::Float64 # minimal difference between current best fitness and second-best one
  fitness_tol::Float64  # minimal difference between current best fitness and the known optmimum

  num_steps::Int       # optimization steps done
  num_better::Int      # number of steps that improved best fitness
  num_better_since_last_report::Int # ditto, since last call to trace_progress()
  last_num_fevals::Int # the number of function evals on the previous step
  num_steps_without_fevals::Int # the number of steps without the function evals

  start_time::Float64 # time optimization started, 0 if not running yet
  stop_time::Float64  # time optimization stopped, 0 if still running
  last_report_time::Float64 # last time trace_progress() was called

  stop_reason::ASCIIString # the reason for algorithm termination, empty if it's not terminated
end

function OptController{O<:Optimizer, E<:Evaluator}(optimizer::O, evaluator::E, params)
  OptController{O,E}(optimizer, evaluator,
        [params[key] for key in Symbol[:ShowTrace, :SaveTrace, :TraceInterval,
                      :MaxSteps, :MaxFuncEvals, :MaxNumStepsWithoutFuncEvals, :MaxTime,
                      :MinDeltaFitnessTolerance, :FitnessTolerance]]...,
        0, 0, 0, 0, 0, 0.0, 0.0, 0.0, "")
end

# stepping optimizer has it's own evaluator, get a reference
OptController(optimizer::SteppingOptimizer, problem::OptimizationProblem, params) = OptController(optimizer, evaluator(optimizer), params)
# all other optimizers are using ProblemEvaluator by default
OptController(optimizer::Optimizer, problem::OptimizationProblem, params) = OptController(optimizer, ProblemEvaluator(problem), params)

# logging/tracing
function tr(ctrl::OptController, msg::AbstractString, obj = nothing)
  if ctrl.show_trace
    print(msg)
    if obj != nothing
      showcompact(obj)
    end
  end
  if ctrl.save_trace
    # No saving for now
  end
end

evaluator(ctrl::OptController) = ctrl.evaluator
problem(ctrl::OptController) = problem(evaluator(ctrl))
isstarted(ctrl::OptController) = ctrl.start_time > 0
isrunning(ctrl::OptController) = isstarted(ctrl) && ctrl.stop_time == 0
isstopped(ctrl::OptController) = ctrl.stop_time > 0

num_steps(ctrl::OptController) = ctrl.num_steps
num_func_evals(ctrl::OptController) = num_evals(ctrl.evaluator)
stop_reason(ctrl::OptController) = ctrl.stop_reason

function elapsed_time(ctrl::OptController)
  isrunning(ctrl) ? time() - ctrl.start_time : ctrl.stop_time - ctrl.start_time
end

function check_stop_condition(ctrl::OptController)
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

    if delta_fitness(ctrl.evaluator.archive) < ctrl.min_delta_fitness_tol
      return "Delta fitness ($(delta_fitness(ctrl.evaluator.archive)) below tolerance ($(ctrl.min_delta_fitness_tol))"
    end

    if fitness_is_within_ftol(ctrl.evaluator, ctrl.fitness_tol)
      return "Fitness ($(best_fitness(ctrl)) within tolerance ($(ctrl.fitness_tol)) of optimum"
    end

    return "" # empty string, no termination
end

function trace_progress(ctrl::OptController)
  ctrl.num_better += ctrl.num_better_since_last_report
  ctrl.last_report_time = time()

  # Always print step number, num fevals and elapsed time
  tr(ctrl, @sprintf("%.2f secs, %d evals, %d steps",
      elapsed_time(ctrl), num_func_evals(ctrl), num_steps(ctrl)))

  # Only print if this optimizer reports on number of better. They return 0
  # if they do not.
  if ctrl.num_better_since_last_report > 0
    tr(ctrl, @sprintf(", improv/step: %.3f (last = %.4f)",
        ctrl.num_better/num_steps(ctrl),
        ctrl.num_better_since_last_report/num_steps(ctrl)))
    ctrl.num_better_since_last_report = 0
  end

  # Always print fitness if num_evals > 0
  if num_func_evals(ctrl) > 0
    tr(ctrl, @sprintf(", %.9f", best_fitness(ctrl)))
  end

  tr(ctrl, "\n")
end

# The ask and tell interface is more general since you can mix and max
# elements from several optimizers using it. However, in this top-level
# execution function we do not make use of this flexibility...
function step!{O<:AskTellOptimizer}(ctrl::OptController{O})
  candidates = ask(ctrl.optimizer)
  rank_by_fitness!(ctrl.evaluator, candidates)
  return tell!(ctrl.optimizer, candidates)
end

# step for SteppingOptimizers
function step!{O<:SteppingOptimizer}(ctrl::OptController{O})
  step!(ctrl.optimizer)
  return 0
end

function run!(ctrl::OptController)
  tr(ctrl, "Starting optimization with optimizer $(name(ctrl.optimizer))\n")
  # FIXME all counters need to be reset if OptController could be run 2nd time
  ctrl.start_time = time()
  ctrl.num_steps = 0
  while isempty(ctrl.stop_reason)
    # Report on progress every now and then...
    if (time() - ctrl.last_report_time) > ctrl.trace_interval
      trace_progress(ctrl)
    end

    # update the counters
    ctrl.num_better_since_last_report += step!(ctrl)
    ctrl.num_steps += 1
    if num_evals(ctrl.evaluator) == ctrl.last_num_fevals
        ctrl.num_steps_without_fevals += 1
    else
        ctrl.num_steps_without_fevals = 0
    end
    ctrl.last_num_fevals = num_func_evals(ctrl)

    ctrl.stop_reason = check_stop_condition(ctrl)
  end
  ctrl.stop_time = time()
  ctrl.num_better += ctrl.num_better_since_last_report
  ctrl.num_better_since_last_report = 0
  tr(ctrl, "\nOptimization stopped after $(ctrl.num_steps) steps and $(elapsed_time(ctrl)) seconds\n")
end

best_candidate(ctrl::OptController) = best_candidate(ctrl.evaluator.archive)
best_fitness(ctrl::OptController) =  best_fitness(ctrl.evaluator.archive)

function show_report(ctrl::OptController, population_stats=false)
  final_elapsed_time = elapsed_time(ctrl)
  tr(ctrl, "Termination reason: $(ctrl.stop_reason)\n")
  tr(ctrl, "Steps per second = $(num_steps(ctrl)/final_elapsed_time)\n")
  tr(ctrl, "Function evals per second = $(num_func_evals(ctrl)/final_elapsed_time)\n")
  tr(ctrl, "Improvements/step = $(ctrl.num_better/ctrl.max_steps)\n")
  tr(ctrl, "Total function evaluations = $(num_func_evals(ctrl))\n")

  if population_stats && isa(ctrl.optimizer, PopulationOptimizer)
    tr(ctrl, "\nMean value (in population) per position:",  mean(population(ctrl.optimizer),1))
    tr(ctrl, "\n\nStd dev (in population) per position:", std(population(ctrl.optimizer),1))
  end

  tr(ctrl, "\n\nBest candidate found: ", best_candidate(ctrl))
  tr(ctrl, "\n\nFitness: ", best_fitness(ctrl))
  tr(ctrl, "\n\n")
end

function write_result(ctrl::OptController, filename = "")
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
