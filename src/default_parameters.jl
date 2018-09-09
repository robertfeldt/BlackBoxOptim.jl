"""
Default parameters for all convenience methods that are exported to the end user.

See `OptRunController` for the description.
"""
const DefaultParameters = ParamsDict(
    :NumDimensions  => :NotSpecified, # Dimension of problem to be optimized
    :SearchRange    => (-1.0, 1.0), # Default search range to use per dimension unless specified
    :SearchSpace    => false, # Search space can be directly specified and will then take precedence over NumDimensions and SearchRange.
    :FitnessScheme  => MinimizingFitnessScheme, # fitness scheme to be used
    :TargetFitness => nothing, # optimal (target) fitness, if known

    :Method => :adaptive_de_rand_1_bin_radiuslimited,

    :MaxTime        => 0.0,
    :MaxFuncEvals   => 0,
    :MaxSteps       => 10000,
    :MaxStepsWithoutProgress => 10000,
    :MinDeltaFitnessTolerance => 1e-50,
    :FitnessTolerance => 1e-8,

    :MaxNumStepsWithoutFuncEvals => 100,

    :NumRepetitions => 1,     # Number of repetitions to run for each optimizer for each problem

    :TraceMode      => :compact,  # Print tracing information during the optimization
    :TraceInterval  => 0.50,  # Minimum number of seconds between consecutive trace messages printed to STDOUT
    :SaveTrace      => false,
    :SaveFitnessTraceToCsv => false, # Save a csv file with information about the major fitness improvement events (only the first event in each fitness magnitude class is saved)
    :SaveParameters => false, # Save parameters to a json file for later scrutiny

    :CallbackFunction => x -> x, # Function to callback to, here just the identity function.
    :CallbackInterval  => 0.0,  # Minimum number of seconds between consecutive callbacks. If 0.0 we never callback.

    :RandomizeRngSeed => true, # Randomize the RngSeed value before using any random numbers.
    :RngSeed        => 1234,   # The specific random seed to set before any random numbers are generated. The seed is randomly selected if RandomizeRngSeed is true, and this parameter is updated with its actual value.

    :PopulationSize => 50
)

function check_and_create_search_space(params::Parameters)
    # If an explicit SearchSpace has been given we use that. It takes precedence.
    if haskey(params, :SearchSpace)
        ss = params[:SearchSpace]
        if isa(ss, SearchSpace)
            return ss
        elseif isa(ss, typeof([(0.0, 1.0)]))
            return RangePerDimSearchSpace(ss)
        elseif ss == false
            # silently fallthrough to the other means of search space specification
        else
            throw(ArgumentError("Using $(typeof(ss)) for SearchSpace is not supported"))
        end
    end

    # If no explicit search space has been given we must create one from the
    # other related parameters.
    if haskey(params, :SearchRange)
        sr = params[:SearchRange]
        # Check that a valid search range has been stated and create the search_space
        # based on it, or bail out.
        if isa(sr, typeof((0.0, 1.0)))
            ndim = params[:NumDimensions]
            if ndim == :NotSpecified
                throw(ArgumentError("You MUST specify NumDimensions= in a solution when giving a SearchRange=$(sr)"))
            end
            return symmetric_search_space(params[:NumDimensions], sr)
        elseif isa(sr, typeof([(0.0, 1.0)]))
            return RangePerDimSearchSpace(sr)
        else
            throw(ArgumentError("Using $(typeof(sr)) for SearchRange is not supported."))
        end
    end

    # No valid search space have been specified => bail.
    throw(ArgumentError("Invalid search space specification ("*
                        "SearchSpace = $(params[:SearchSpace]), "*
                        "SearchRange = $(params[:SearchRange]), "*
                        "NumDimensions = $(params[:NumDimensions]))"))
end

function check_valid!(params::Parameters)
    # Check that max_time is larger than zero if it has been specified.
    if haskey(params, :MaxTime)
        if !isa(params[:MaxTime], Number) || params[:MaxTime] < 0.0
            throw(ArgumentError("MaxTime parameter must be a non-negative number"))
        elseif params[:MaxTime] > 0.0
            params[:MaxTime] = convert(Float64, params[:MaxTime])
            params[:MaxFuncEvals] = 0
            params[:MaxSteps] = 0
        end
    end

    # Check that a valid number of fevals has been specified. Print warning if higher than 1e8.
    if haskey(params,:MaxFuncEvals)
        if !isa(params[:MaxFuncEvals], Integer) || params[:MaxFuncEvals] < 0.0
            throw(ArgumentError("MaxFuncEvals parameter MUST be a non-negative integer"))
        elseif params[:MaxFuncEvals] > 0.0
            if params[:MaxFuncEvals] >= 1e8
                warn("Number of allowed function evals is $(params[:MaxFuncEvals]); this can take a LONG   time")
            end
            params[:MaxFuncEvals] = convert(Int, params[:MaxFuncEvals])
            params[:MaxSteps] = 0
        end
    end

    # Check that a valid number of iterations has been specified. Print warning if higher than 1e8.
    if haskey(params, :MaxSteps)
        if !isa(params[:MaxSteps], Number) || params[:MaxSteps] < 0.0
            throw(ArgumentError("The number of iterations (MaxSteps) MUST be a non-negative number"))
        elseif params[:MaxSteps] > 0.0
            if params[:MaxSteps] >= 1e8
                warn("Number of allowed iterations is $(params[:MaxSteps]); this can take a LONG time")
            end
            params[:MaxSteps] = convert(Int, params[:MaxSteps])
        end
    end

    # Check that a valid population size has been given.
    if params[:PopulationSize] < 2
        # FIXME why? What if we use popsize of 1 for optimizers that improve on one solution?
        throw(ArgumentError("The population size MUST be at least 2"))
    end

    method = params[:Method]
    # Check that a valid method has been specified and then set up the optimizer
    if !isa(method, Symbol) || method ∉ BlackBoxOptim.MethodNames
        throw(ArgumentError("The method specified, $(method), is NOT among the valid methods:   $(MethodNames)"))
    end
end
