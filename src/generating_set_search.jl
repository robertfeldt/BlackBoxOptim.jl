# Generating Set Search as described in Kolda2003:
#  Kolda, Tamara G., Robert Michael Lewis, and Virginia Torczon. "Optimization
#  by direct search: New perspectives on some classical and modern methods."
#  SIAM review 45.3 (2003): 385-482.
#

"""
A supertype for all generating set searcher-like algorithms.
"""
abstract type DirectSearcher <: SteppingOptimizer end

"""
`DirectionGenerator` generates the search directions to use at each step of
a GSS search.
"""
abstract type DirectionGenerator end

struct ConstantDirectionGen <: DirectionGenerator
    directions::Matrix{Float64}

    ConstantDirectionGen(directions) = new(directions)
end

directions_for_k(cg::ConstantDirectionGen, k) =
    cg.directions # For a ConstantDirectionGen it is always the same regardless of k...

# We can easily do a compass search with GSS by generating directions
# individually (+ and -) for each coordinate.
compass_search_directions(n) = ConstantDirectionGen([Matrix{Float64}(I, n,n) -Matrix{Float64}(I, n, n)])

const GSSDefaultParameters = ParamsDict(
    :DeltaTolerance => 1e-10,       # GSS has converged if the StepSize drops below this tolerance level
    :InitialStepSizeFactor => 0.50,        # Factor times the minimum search space diameter to give the initial StepSize
    :RandomDirectionOrder => true,  # Randomly shuffle the order in which the directions are used for each step
    :StepSizeGamma => 2.0,          # Factor by which step size is multiplied if improved point is found. Should be >= 1.0.
    :StepSizePhi => 0.5,            # Factor by which step size is multiplied if NO improved point is found. Should be < 1.0.
    :StepSizeMax => prevfloat(typemax(Float64)) # A limit on the step size can be set but is typically not => Inf.
)

calc_initial_step_size(ss, stepSizeFactor = 0.5) = stepSizeFactor * minimum(diameters(ss))

"""
Generating Set Search as described in Kolda2003:
  Kolda, Tamara G., Robert Michael Lewis, and Virginia Torczon. "Optimization
  by direct search: New perspectives on some classical and modern methods."
  SIAM review 45.3 (2003): 385-482.
"""
mutable struct GeneratingSetSearcher{V<:Evaluator, D<:DirectionGenerator, E<:EmbeddingOperator} <: DirectSearcher
    direction_gen::D
    evaluator::V
    embed::E
    search_space::SearchSpace
    n::Int
    k::Int
    random_dir_order::Bool       # shuffle the order of directions?
    step_size_factor::Float64    # initial step size factor
    step_size_gamma::Float64     # step size factor if improved
    step_size_phi::Float64       # step size factor if no improvement
    step_size_max::Float64       # maximal step size
    step_tol::Float64            # step delta tolerance
    step_size::Float64           # current step size
    x::Individual
    xfitness::Float64

    function GeneratingSetSearcher(
            evaluator::V, dgen::D, embed::E,
            random_dir_order::Bool, step_size_factor::Float64, step_size_gamma::Float64,
            step_size_phi::Float64, step_size_max::Float64,
            step_tol::Float64) where {V<:Evaluator, D<:DirectionGenerator, E<:EmbeddingOperator}
        # FIXME check parameters ranges
        n = numdims(evaluator)
        ss = search_space(evaluator)
        x = rand_individual(ss)
        new{V, D, E}(dgen, evaluator, embed, ss, n, 0,
            random_dir_order, step_size_factor,
            step_size_gamma, step_size_phi,
            step_size_max, step_tol,
            calc_initial_step_size(ss, step_size_factor),
            x, fitness(x, evaluator))
    end
end

# by default use RandomBound embedder
function GeneratingSetSearcher(problem::OptimizationProblem, parameters::Parameters)
    params = chain(GSSDefaultParameters, parameters)
    GeneratingSetSearcher(ProblemEvaluator(problem),
                        get(params, :DirectionGenerator, compass_search_directions(numdims(problem))),
                        RandomBound(search_space(problem)),
                        params[:RandomDirectionOrder],
                        params[:InitialStepSizeFactor], params[:StepSizeGamma],
                        params[:StepSizePhi], params[:StepSizeMax],
                        params[:DeltaTolerance])
end

# We also include the name of the direction generator.
name(opt::GeneratingSetSearcher) = "GeneratingSetSearcher($(typeof(opt.direction_gen)))"

has_converged(gss::GeneratingSetSearcher) = gss.step_size < gss.step_tol

function step!(gss::GeneratingSetSearcher)
    if has_converged(gss)
        # Restart from a random point
        gss.x = rand_individual(gss.search_space)
        gss.xfitness = fitness(gss.x, gss.evaluator)
        gss.step_size = calc_initial_step_size(gss.search_space, gss.step_size_factor)
    end

    # Get the directions for this iteration
    gss.k += 1
    directions = directions_for_k(gss.direction_gen, gss.k)

    # Set up order vector from which we will take the directions after possibly shuffling it
    order = collect(1:size(directions, 2))
    if gss.random_dir_order
        shuffle!(order)
    end

    # Check all directions to find a better point; default is that no one is found.
    found_better = false
    candidate = zeros(gss.n, 1)

    # Loop over directions until we find an improvement (or there are no more directions to check).
    for direction in order
        candidate = gss.x + gss.step_size .* directions[:, direction]
        apply!(gss.embed, candidate, gss.x)

        if is_better(candidate, gss.xfitness, gss.evaluator)
            found_better = true
            break
        end
    end

    if found_better
        gss.x = candidate
        gss.xfitness = last_fitness(gss.evaluator)
        gss.step_size *= gss.step_size_gamma
    else
        gss.step_size *= gss.step_size_phi
    end
    gss.step_size = min(gss.step_size, gss.step_size_max)

    return gss
end
