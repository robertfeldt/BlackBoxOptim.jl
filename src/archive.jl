"""
`Archive` saves information about promising solutions during an optimization run.
"""
abstract type Archive{F,FS<:FitnessScheme} end

numdims(a::Archive) = a.numdims
fitness_type(a::Archive{F}) where F = F
fitness_scheme(a::Archive) = a.fit_scheme

"""
    archived_fitness(fit, a::Archive)

Converts given fitness `fit` to the format used by the archive `a`.
"""
archived_fitness(fit::Any, a::Archive) = convert(fitness_type(a), fit, fitness_scheme(a))

"""
Base class for individuals stored in `Archive`.
"""
abstract type ArchivedIndividual{F} <: FitIndividual{F} end

tag(indi::ArchivedIndividual) = indi.tag

"""
Individual stored in `TopListArchive`.
"""
struct TopListIndividual{F} <: ArchivedIndividual{F}
    params::Individual
    fitness::F
    tag::Int

    TopListIndividual(params::AbstractIndividual, fitness::F, tag::Int) where F =
        new{F}(params, fitness, tag)
end

Base.:(==)(x::TopListIndividual{F}, y::TopListIndividual{F}) where F =
  (x.fitness == y.fitness) && (x.params == y.params)

"""
Fitness as stored in `TopListArchive`.
"""
struct TopListFitness{F<:Number}
    fitness::F            # current fitness
    fitness_improvement_ratio::Float64
    num_fevals::Int     # number of fitness evaluations so far
    timestamp::Float64    # when archived
end

fitness(f::TopListFitness) = f.fitness

"""
Archive that maintains a top list of the best performing (best fitness)
candidates seen so far.


The two last best fitness values could be used to estimate a confidence interval for how much improvement
potential there is.
"""
mutable struct TopListArchive{F<:Number,FS<:FitnessScheme} <: Archive{F,FS}
    fit_scheme::FS        # Fitness scheme used
    start_time::Float64   # Time when archive created, we use this to approximate the starting time for the opt...
    numdims::Int          # Number of dimensions in the optimization problem. Needed for confidence interval estimation.

    num_fitnesses::Int    # Number of calls to add_candidate

    capacity::Int         # Max size of top lists
    candidates::Vector{TopListIndividual{F}}  # Top candidates and their fitness values

    # Stores a fitness history that we can later print to a csv file.
    # For each magnitude class (as defined by `magnitude_class()` function below) we
    # we save the first entry of that class. The tuple saved for each magnitude
    # class is: `(magnitude_class, time, num_fevals, fitness, width_of_confidence_interval)`
    fitness_history::Vector{TopListFitness{F}}

    function TopListArchive(fit_scheme::FS, numdims::Integer,
                            capacity::Integer = 10) where FS<:FitnessScheme
        F = fitness_type(FS)
        new{F,FS}(fit_scheme, time(), numdims, 0, capacity,
                  TopListIndividual{F}[], TopListFitness{F}[])
    end
end

fitness_scheme(a::TopListArchive) = a.fit_scheme
capacity(a::TopListArchive) = a.capacity
Base.length(a::TopListArchive) = length(a.candidates)
Base.isempty(a::TopListArchive) = isempty(a.candidates)

best_candidate(a::TopListArchive) = a.candidates[1].params
best_fitness(a::TopListArchive) = !isempty(a.candidates) ? fitness(a.candidates[1]) : nafitness(fitness_scheme(a))
last_top_fitness(a::TopListArchive) = !isempty(a.candidates) ? fitness(a.candidates[end]) : nafitness(fitness_scheme(a))

"""
    delta_fitness(a::TopListArchive)

The difference between the current best fitness and the former best fitness.
"""
function delta_fitness(a::TopListArchive)
    if length(a.fitness_history) < 2
        Inf
    else
        # FIXME aggregate fitness?
        abs(a.fitness_history[end].fitness - a.fitness_history[end-1].fitness)
    end
end

function check_stop_condition(a::TopListArchive, p::OptimizationProblem, ctrl)
    if delta_fitness(a) < ctrl.min_delta_fitness_tol
        return "Delta fitness ($(delta_fitness(a))) below tolerance ($(ctrl.min_delta_fitness_tol))"
    end

    if fitness_is_within_ftol(p, best_fitness(a), ctrl.fitness_tol)
        return "Fitness ($(best_fitness(ctrl))) within tolerance ($(ctrl.fitness_tol)) of optimum"
    end

    return "" # no conditions met
end

"""
    add_candidate!(a::TopListArchive, fitness::F, candidate[, tag=0][, num_fevals=-1])

Add a candidate with a fitness to the archive (if it is good enough).
"""
function add_candidate!(a::TopListArchive{F}, fitness::F, candidate::AbstractIndividual,
                        tag::Int=0, num_fevals::Int=-1) where F
    a.num_fitnesses += 1
    if (num_fevals == -1) num_fevals = a.num_fitnesses end

    if isempty(a.fitness_history) || is_better(fitness, best_fitness(a), fitness_scheme(a))
        # Save fitness history so we can reconstruct the most important events later.
        push!(a.fitness_history, TopListFitness{F}(fitness, fitness_improvement_ratio(a, fitness), num_fevals, time()))
    end

    if length(a) < capacity(a) ||
       !isempty(a.candidates) && is_better(fitness, last_top_fitness(a), fitness_scheme(a))
        new_cand = TopListIndividual(copy(candidate), fitness, tag)
        fs = fitness_scheme(a)
        ix = searchsortedfirst(a.candidates, new_cand,
                               # FIXME use lt=fitness_scheme(a) when v0.5 #14919 would be fixed
                               by=BlackBoxOptim.fitness, lt=(x, y) -> is_better(x, y, fs))
        if ix > length(a) || a.candidates[ix] != new_cand
            insert!(a.candidates, ix, new_cand)
        end
        (length(a) > capacity(a)) && pop!(a.candidates) # don't grow over the capacity
    end
end

"""
Get a tuple of the sign and the magnitude of the value rounded to the first digit.
Used for archiving candidates separately for each magnitude class.
"""
function magnitude_class(f)
    f = float(f)
    if f == 0.0
        (-1.0, 1e100)
    else
        (sign(f), floor(10.0*log10(abs(f)))/10.0)
    end
end

function fitness_improvement_ratio(a::Archive, newFitness)
    try
        bestfitness = best_fitness(a)
        return abs( (bestfitness - newFitness) / bestfitness )
    catch
        return 0.0
    end
end

"""
Get the distance from a fitness value to the optimum/best known fitness value.
"""
distance_to_optimum(fitness, bestfitness) = abs(fitness - bestfitness)

function fitness_history_csv_header(a::Archive)
    "Date,Time,ElapsedTime,Magnitude,NumFuncEvals,FitnessImprovementRatio,Fitness"
end

function save_fitness_history_to_csv_file(a::Archive, filename = "fitness_history.csv";
    header_prefix = "", line_prefix = "",
    include_header = true, bestfitness = nothing)

    fh = open(filename, "a+")

    if include_header
        header = [header_prefix, fitness_history_csv_header(a)]

        if bestfitness != nothing
        push!(header, "DistanceToBestFitness")
        end

        println(fh, join(header, ","))
    end

    for af in a.fitness_history
        mc = magnitude_class(af.fitness)

        line = [line_prefix, strftime("%Y-%m-%d,%T", af.timestamp),
            af.timestamp-a.start_time,
            mc[1]*mc[2], af.num_fevals, af.fitness_improvement_ratio, af.fitness]

        if bestfitness != nothing
            push!(line, distance_to_optimum(af.fitness, bestfitness))
        end

        println(fh, join(line, ","))
    end

    close(fh)
end

"""
    merge_fitness_histories(histories)

Merge the collection of multiple fitness histories and calculate the `min`, `max`, `avg` and `median`
values for fitness and FIR (fitness improvement ratio) at each point where the fitness is changing.
"""
function merge_fitness_histories(histories)
    numhist = length(histories)
    counts = ones(numhist)
    fitnesses = zeros(numhist)
    firs = zeros(numhist)
    current_feval = 1

    # Find max history length
    maxlen = mapreduce( length, max, histories )

    while maximum(counts) < maxlen
        # Find min feval for current events, this is the next current feval
        for i in 1:numhist
            t, nf, f, fir = histories[i][counts[i]]
            # FIXME ???
        end
    end
end

# Merge the fitness histories and save the average values of the fitness,
# and distance to best fitness for each change in any of the histories.
#function merge_fitness_histories_to_csv(archives::Archive[], filename = "fitness_history.csv";
#  header_prefix = "", line_prefix = "",
#  include_header = true, bestfitness = nothing)
#
#end

"""
Calculate the width of the confidence interval at a certain `p`-value.
This is based on the paper:
    Carvalho (2011), "Confidence intervals for the minimum of a
    function using extreme value statistics"

This means that the current estimate of the confidence interval for the minimum
of the optimized function lies within the interval

```
] l1 - w, l1 [
```

with probability ``(1-p)`` as the number of sampled function points goes to infinity,
where
```
w = width_of_confidence_interval(a, p)
l1 = best_fitness(a)
```
"""
function width_of_confidence_interval(a::Archive, p = 0.01)
    if length(a) < 2
        return Inf
    else
        l1 = fitness(a.candidates[1])
        l2 = fitness(a.candidates[2])
        # We use abs below so it works also for maximization.
        abs(l2 - l1) / ( (1-p)^(-2/a.numdims) - 1)
    end
end

"""
    fitness_improvement_potential(a::Archive[, p = 0.01])

Calculate the solution improvement potential.

The percentage improvement that can be expected from the current fitness value at a given `p`-value.
In theory, an optimization run should be terminated when this value is very
small, i.e. there is little improvement potential left in the run.
"""
fitness_improvement_potential(a::Archive, p = 0.01) =
    width_of_confidence_interval(a, p) / abs(best_fitness(a))
