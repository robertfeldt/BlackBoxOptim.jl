"""
    setup_problem(problem, parameters::Parameters)

Set up a fixed-dimensional optimization problem.

There are several `setup_problem()` method that accept different type of
`problem` argument:
    * `OptimizationProblem`
    * function (`:NumDimensions` has to be specified in `parameters`)
    * `FunctionBasedProblemFamily` (`:NumDimensions` has to be specified in `parameters`)
"""
function setup_problem(problem::OptimizationProblem, parameters::Parameters)
    return problem, parameters
end

function setup_problem(family::FunctionBasedProblemFamily, parameters::Parameters)
    # If an anydim problem was given the dimension param must have been specified.
    if parameters[:NumDimensions] == :NotSpecified
        throw(ArgumentError("You MUST specify NumDimensions= when a problem family is given"))
    end
    problem = instantiate(family, parameters[:NumDimensions])

    return problem, parameters
end

function setup_problem(func::Function, parameters::Parameters)
    ss = check_and_create_search_space(parameters)

    # Now create an optimization problem with the given information. We currently reuse the type
    # from our pre-defined problems so some of the data for the constructor is dummy.
    problem = FunctionBasedProblem(func, "<unknown>", parameters[:FitnessScheme], ss,
                                   parameters[:TargetFitness])

    # validate fitness: create a random solution from the search space and ensure that fitness(problem) returns fitness_type(problem).
    ind = rand_individual(search_space(problem))
    res = fitness(ind, problem)
    if !isa(res, fitness_type(problem))
      throw(ArgumentError("The supplied fitness function does NOT return the expected fitness type "*
                          "when called with a potential solution "*
                          "(when called with $(ind) it returned $(res)) so we cannot optimize it!"))
    end

    return problem, parameters
end

"""
    bboptimize(problem[; parameters::Associative, kwargs...])

Solve given optimization `problem`.

See `setup_problem()` for the types of `problem` supported.
In addition, the `problem` could be `OptController` containing the
results of the previous optimization runs.

The optimization method parameters could be specified via `kwargs` or `parameters` arguments.

Returns `OptimizationResults` instance.

See also `bbsetup()`.
"""
function bboptimize(optctrl::OptController; kwargs...)
    if length(kwargs) > 0
        update_parameters!(optctrl, kwargs2dict(kwargs))
    end
    run!(optctrl)
end

function bboptimize(functionOrProblem, parameters::Parameters = EMPTY_PARAMS; kwargs...)
    optctrl = bbsetup(functionOrProblem, parameters; kwargs...)
    run!(optctrl)
end

"""
    bbsetup(problem[; parameters::Associative, kwargs...])

Set up optimization method for a given problem.

See `setup_problem()` for the types of `problem` supported.
The optimization method parameters could be specified via `kwargs` or `parameters` arguments.

Returns the initialized `OptController` instance. To actually run the method
call `bboptimize()` or `run!()`.
"""
function bbsetup(functionOrProblem, parameters::Parameters = EMPTY_PARAMS; kwargs...)
    parameters = chain(convert(ParamsDict, parameters), kwargs2dict(kwargs))
    problem, params = setup_problem(functionOrProblem, chain(DefaultParameters, parameters))
    check_valid!(params)

    optimizer_func = chain(SingleObjectiveMethods, MultiObjectiveMethods)[params[:Method]]
    optimizer = optimizer_func(problem, params)

    # Now set up a controller for this problem. This will handle
    # application of optimizer, checking for termination conditions
    # as well as collecting the statistics about the number of function evals,
    # keep an archive and top list of candidates.
    ctrl = OptController(optimizer, problem, params)

    return ctrl
end
