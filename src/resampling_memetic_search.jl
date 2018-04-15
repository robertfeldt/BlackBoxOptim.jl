const RSDefaultParameters = ParamsDict(
    :PrecisionRatio    => 0.40, # 40% of the diameter is used as the initial step length
    :PrecisionTreshold => 1e-6  # They use 1e-6 in the paper.
)

"""
The variants of the memetic search algorithms RS and RIS.
However, we have modified them since they did not give very good performance when
implemented as described in the papers below. Possibly, the papers are not
unambigous and I have misinterpreted something from them...

The "Resampling Search" (RS) memetic algorithm is described in:

    F. Caraffini, F. Neri, M. Gongora and B. N. Passow, "Re-sampling Search: A
    Seriously Simple Memetic Approach with a High Performance", 2013.

and its close sibling "Resampling Inheritance Search" (RIS) is described in:

    F. Caraffini, F. Neri, B. N. Passow and G. Iacca, "Re-sampled Inheritance
    Search: High Performance Despite the Simplicity", 2013.
"""
mutable struct ResamplingMemeticSearcher{E<:Evaluator} <: SteppingOptimizer
    name::String
    params::Parameters
    evaluator::E
    resampling_func::Function

    precisions      # Cache the starting precision values so we need not calc them for each step
    diameters       # Cache the diameters...

    elite           # Current elite (best) candidate
    elite_fitness   # Fitness of current elite

    # Constructor for RS:
    function ResamplingMemeticSearcher(evaluator::E, parameters::Parameters,
            resampling_function::Function, name::String) where {E<:Evaluator}
        params = chain(RSDefaultParameters, parameters)
        elite = rand_individual(search_space(evaluator))
        diams = diameters(search_space(evaluator))
        new{E}(name, params, evaluator, resampling_function,
               params[:PrecisionRatio] * diams, diams,
               elite, fitness(elite, evaluator))
    end
end

name(rs::ResamplingMemeticSearcher) = rs.name

const RISDefaultParameters = ParamsDict(
    :InheritanceRatio => 0.30   # On average, 30% of positions are inherited when resampling in RIS
)

ResamplingMemeticSearcher(problem::OptimizationProblem,
        params::Parameters = EMPTY_PARAMS,
        resampling_function = random_resample,
        name = "Resampling Memetic Search (RS)") =
    ResamplingMemeticSearcher(ProblemEvaluator(problem), params, resampling_function, name)

# Constructor for the RIS:
ResamplingInheritanceMemeticSearcher(problem::OptimizationProblem, parameters::Parameters = EMPTY_PARAMS) =
    ResamplingMemeticSearcher(problem,
            chain(RSDefaultParameters, RISDefaultParameters, parameters),
            random_resample_with_inheritance,
            "Resampling Inheritance Memetic Search (RIS)")

resampling_memetic_searcher(problem::OptimizationProblem, params::Parameters) =
    ResamplingMemeticSearcher(problem, params)

resampling_inheritance_memetic_searcher(problem::OptimizationProblem, params::Parameters) =
    ResamplingInheritanceMemeticSearcher(problem, params)

# For Resampling Search (RS) the resample is purely random.
random_resample(rms::ResamplingMemeticSearcher) = rand_individual(search_space(rms.evaluator))

# For Resampling Inheritance Search (RIS) the resample has an inheritance component.
function random_resample_with_inheritance(rms::ResamplingMemeticSearcher)
    xt = random_resample(rms)
    n = numdims(rms.evaluator)
    i = rand(1:n)
    Cr = 0.5^(1/(rms.params[:InheritanceRatio]*n)) # See equation 3 in the RIS paper
    k = 1

    while rand() <= Cr && k < n
        xt[i] = rms.elite[i]
        i = 1 + mod(i, n)
        k += 1
    end

    return xt
end

function step!(rms::ResamplingMemeticSearcher)
    # First randomly sample two candidates and select the best one. It seems
    # RS and RIS might be doing this in two different ways but use the RS way for
    # now.
    trial, fitness = best_of(rms.resampling_func(rms), rms.resampling_func(rms), rms.evaluator)

    # Update elite if new trial is better. This is how they write it in the RIS paper
    # but in the RS paper it seems they always update the elite. Unclear! To me it
    # seems we should always update since we have already done a local search from
    # the current elite so it has a "head start" compared to new sampled points
    # which have not yet gone through local refinement. Since the evaluator/archive
    # keeps the best candidates anyway there is no risk for us in always overwriting the elite...
    set_as_elite_if_better(rms, trial, fitness)

    # Then run the local search on the elite one until step length too small.
    trial, fitness = local_search(rms, trial, fitness)
    set_as_elite_if_better(rms, trial, fitness)

    return trial, fitness
end

function set_as_elite_if_better(rms::ResamplingMemeticSearcher, candidate, fitness)
    if is_better(fitness, rms.elite_fitness, rms.evaluator)
        rms.elite = candidate
        rms.elite_fitness = fitness
        return true
    else
        return false
    end
end

stop_due_to_low_precision(rms::ResamplingMemeticSearcher, precisions) =
    norm(precisions ./ rms.diameters) < rms.params[:PrecisionTreshold]

function local_search(rms::ResamplingMemeticSearcher, xstart, fitness)
    ps = copy(rms.precisions)
    xt = copy(xstart)
    tfitness = copy(fitness)

    searchSpace = search_space(rms.evaluator)
    ssmins, ssmaxs = mins(searchSpace), maxs(searchSpace)
    n = numdims(rms.evaluator)

    indices = collect(1:n)

    while !stop_due_to_low_precision(rms, ps)
        xs = copy(xt)
        # We randomize the order that each decision var is changed. This is not done in the orig papers.
        shuffle!(indices)

        for j in 1:n
            i = indices[j]

            # This is how it is written in orig papers. To me it seems better to
            # take the step in a random direction; why prioritize one direction?
            xs[i] = xt[i] - ps[i]

            # rand bound from target if below min
            if xs[i] < ssmins[i]
                xs[i] = ssmins[i] + rand() * (xt[i] - ssmins[i])
            end

            if is_better(xs, tfitness, rms.evaluator)
                xt[i] = xs[i]
                tfitness = last_fitness(rms.evaluator)
            else
                xs[i] = xt[i] + ps[i]/2

                # rand bound from target if above max
                if xs[i] > ssmaxs[i]
                    xs[i] = xt[i] + rand() * (ssmaxs[i] - xt[i])
                end

                if is_better(xs, tfitness, rms.evaluator)
                    xt[i] = xs[i]
                    tfitness = last_fitness(rms.evaluator)
                end
            end
        end

        if is_better(tfitness, fitness, rms.evaluator)
            fitness = tfitness
            xstart = xt
        else
            ps = ps / 2
        end
    end

    return xstart, fitness
end
