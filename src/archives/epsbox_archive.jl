"""
Individual representing the solution from the Pareto set.
"""
struct FrontierIndividual{F} <: ArchivedIndividual{F}
    fitness::F
    params::Individual
    tag::Int                            # tag of the individual (e.g. gen.op. ID)
    num_fevals::Int                     # number of fitness evaluations so far
    n_restarts::Int                     # the number of method restarts so far
    timestamp::Float64                  # when archived

    FrontierIndividual{F}(fitness::F,
                   params, tag, num_fevals, n_restarts, timestamp=time()) where F =
        new{F}(fitness, params, tag, num_fevals, n_restarts, timestamp)

    FrontierIndividual(fitness::F,
                   params, tag, num_fevals, n_restarts, timestamp=time()) where F =
        new{F}(fitness, params, tag, num_fevals, n_restarts, timestamp)
end

tag(indi::FrontierIndividual) = indi.tag

"""
Individual stored in `EpsBoxArchive`.
"""
const EpsBoxFrontierIndividual{N,F<:Number} = FrontierIndividual{IndexedTupleFitness{N,F}}
const FrontierSpatialIndex{N,F<:Number} = SI.SpatialIndex{Int,N,EpsBoxFrontierIndividual{N,F}}
const FrontierRTree{N,F<:Number} = SI.RTree{Int,N,EpsBoxFrontierIndividual{N,F}}

# traits for spatial indexing
SI.mbrtrait(::Type{EpsBoxFrontierIndividual{N,F}}) where {N,F} = SI.HasMBR{SI.Rect{Int,N}}
SI.mbr(indi::EpsBoxFrontierIndividual{N,F}) where {N,F} = SI.Rect(indi.fitness.index, indi.fitness.index)
SI.idtrait(::Type{<:EpsBoxFrontierIndividual}) = SI.HasID{Int}
SI.id(indi::EpsBoxFrontierIndividual) = indi.num_fevals

EpsBoxFrontierIndividual(fitness::IndexedTupleFitness{N,F},
               params, tag, num_fevals, n_restarts, timestamp=time()) where {N,F} =
    FrontierIndividual(fitness, params, tag, num_fevals, n_restarts, timestamp)

"""
ϵ-box archive saves only the solutions that are not ϵ-box
dominated by any other solutions in the archive.

It also counts the number of candidate solutions that have been added
and how many ϵ-box progresses have been made.
"""
mutable struct EpsBoxArchive{N,F,FS<:EpsBoxDominanceFitnessScheme} <: Archive{IndexedTupleFitness{N,F},FS}
    fit_scheme::FS        # Fitness scheme used
    start_time::Float64   # Time when archive created, we use this to approximate the starting time for the opt...

    num_candidates::Int               # Number of calls to add_candidate!()
    best_front_elem::EpsBoxFrontierIndividual{N,F} # best frontier element: the candidate with the best aggregated fitness
    last_progress::Int                # when (wrt num_candidates) last ϵ-progress has occured
    last_restart::Int                 # when (wrt num_dlast) last restart has occured
    n_restarts::Int                   # the counter of the method restarts
    n_oversize_inserts::Int           # how many times the candidates were inserted into oversized archive

    max_size::Int         # maximal frontier size
    # TODO allow different frontier containers?
    # see e.g. Altwaijry & Menai "Data Structures in Multi-Objective Evolutionary Algorithms", 2012
    frontier::FrontierRTree{N,F}  # candidates along the fitness Pareto frontier

    EpsBoxArchive(fit_scheme::EpsBoxDominanceFitnessScheme{N,F};
                  max_size::Integer = 1_000_000,
                  leaf_capacity::Integer = 10,
                  branch_capacity::Integer = 10) where {N,F} =
        new{N,F,typeof(fit_scheme)}(fit_scheme, time(), 0,
                                    EpsBoxFrontierIndividual{N,F}(nafitness(fit_scheme), Individual(), 0, 0, 0, NaN),
                                    0, 0, 0, 0, max_size,
                                    FrontierRTree{N,F}(leaf_capacity=leaf_capacity,
                                                       branch_capacity=branch_capacity))
end

EpsBoxArchive(fit_scheme::EpsBoxDominanceFitnessScheme, params::Parameters) =
    EpsBoxArchive(fit_scheme, max_size=params[:MaxArchiveSize],
                  leaf_capacity=params[:LeafCapacity],
                  branch_capacity=params[:BranchCapacity])

const EpsBoxArchive_DefaultParameters = ParamsDict(
    :MaxArchiveSize => 100_000,
    :LeafCapacity => 10,
    :BranchCapacity => 10
)

Base.eltype(::Type{EpsBoxArchive{N,F}}) where {N,F} = EpsBoxFrontierIndividual{N,F}

Base.length(a::EpsBoxArchive) = length(a.frontier)
Base.isempty(a::EpsBoxArchive) = isempty(a.frontier)
capacity(a::EpsBoxArchive) = a.max_size
numdims(a::EpsBoxArchive) = !isempty(a.frontier) ? length(a.frontier[1].params) : 0

"""
    pareto_frontier(a::EpsBoxArchive)

Get the iterator to the individuals on the Pareto frontier.
"""
pareto_frontier(a::EpsBoxArchive) = a.frontier

"""
    rand_front_elem(a::EpsBoxArchive)

Get random Pareto frontier element.

Returns `nothing` if frontier is empty.
"""
function rand_front_elem(a::EpsBoxArchive)
    isempty(a) && return nothing

    node = a.frontier.root
    while SI.level(node) > 1
        # descend to a random child
        node = node[rand(eachindex(SI.children(node)))]
    end
    return node[rand(eachindex(SI.children(node)))]
end

"""
    noprogress_streak(a::EpsBoxArchive, [since_restart])

Get the number of `add_candidate!()` calls since the last ϵ-progress.
If `since_restart` is specified, the number is relative to the last
restart.
"""
noprogress_streak(a::EpsBoxArchive; since_restart::Bool=false) =
    since_restart ?
        a.num_candidates - max(a.last_progress, a.last_restart) :
        a.num_candidates - a.last_progress

has_best_front_elem(a::EpsBoxArchive) = isfinite(a.best_front_elem.timestamp)
best_front_elem(a::EpsBoxArchive) = has_best_front_elem(a) ? a.best_front_elem : nothing
best_candidate(a::EpsBoxArchive) = has_best_front_elem(a) ? params(a.best_front_elem) : nothing
best_fitness(a::EpsBoxArchive) = has_best_front_elem(a) ? fitness(a.best_front_elem) : nafitness(a.fit_scheme)

function notify!(a::EpsBoxArchive, event::Symbol)
    if event == :restart
        a.n_restarts += 1
        a.last_restart = a.num_candidates
        # rebuild the frontier upon the first restart, or if the frontier has changed significantly
        if a.n_restarts == 1 || a.frontier.nelem_insertions > 3length(a.frontier)
            # bulk-reload R-tree (leave leaves underfilled to postpone reinserts)
            a.frontier = SI.load!(similar(a.frontier), a.frontier,
                                  leaf_fill = floor(Int, a.frontier.fill_factor * SI.capacity(SI.Leaf, a.frontier)))
        end
    end
    return a
end

"""
    tagcounts(a::EpsBoxArchive)

Count the tags of individuals on the ϵ-box frontier.
Each restart the individual remains in the frontier discounts it by `θ`.

Returns the `tag`→`count` dictionary.
"""
function tagcounts(a::EpsBoxArchive, θ::Number = 1.0)
    (0.0 < θ <= 1.0) || throw(ArgumentError("θ ($θ) should be in (0.0, 1.0] range"))
    res = Dict{Int,Float64}()
    for indi in pareto_frontier(a)
        curtag = tag(indi)
        if curtag > 0
            curcounts = get!(res, curtag, 0.0)
            res[curtag] = curcounts+θ^(a.n_restarts-indi.n_restarts)
        end
    end
    return res
end

function add_candidate!(a::EpsBoxArchive{N,F}, cand_fitness::IndexedTupleFitness{N,F},
                        candidate, tag::Int=0, num_fevals::Int=-1) where {N,F}
    a.num_candidates += 1
    if num_fevals == -1
        num_fevals = a.num_candidates
    end
    #@debug "New fitness: $(cand_fitness.orig) agg=$(cand_fitness.agg)"
    #@debug "Params: $candidate"
    # TODO don't traverse the tree twice: some SpatialIndexing API to combine isempty() with the follow-up findfirst() call?
    if !isempty(a.frontier, DominanceCone{is_minimizing(a.fit_scheme)}(cand_fitness.index)) # candidate is dominated by frontier
        if length(a.frontier) <= a.max_size
            a.n_oversize_inserts = 0 # reset the counter since the size is ok
        end
        return a
    end
    # find the Pareto front element with exactly the same indexed fitness as the candidate
    frontix = findfirst(a.frontier, SI.Point(cand_fitness.index))
    if frontix !== nothing # found front element, no eps-progress
        front_leaf, frontel_pos = frontix
        frontel = front_leaf[frontel_pos]
        front_fitness = fitness(frontel)
        hat, index_match = hat_compare(cand_fitness, front_fitness, a.fit_scheme)
        @assert index_match # should have the same indexed fitness since that's how it was found
        if hat < 0 # new fitness dominates (but not eps-dominates) the one in archive, replace the element
            SI.children(front_leaf)[frontel_pos] = frontel =
                EpsBoxFrontierIndividual(cand_fitness, copyto!(frontel.params, candidate), tag, num_fevals, a.n_restarts)
        end
    else # eps-progress: non-dominated solution that has some eps-indices different from all current ones
        a.last_progress = a.num_candidates
        hat = -1
        # remove all dominated frontier elements, if any
        SI.subtract!(a.frontier, DominanceCone{!is_minimizing(a.fit_scheme)}(cand_fitness.index))
        #@debug "Appended non-dominated element to the frontier"
        frontel = EpsBoxFrontierIndividual(cand_fitness, copy(candidate), tag, num_fevals, a.n_restarts)
        insert!(a.frontier, frontel)
        if length(a.frontier) > a.max_size
            a.n_oversize_inserts += 1 # throw(error("Pareto frontier exceeds maximum size"))
        end
    end
    # check if the new candidate has better aggregate score
    if !has_best_front_elem(a) # initialize best element
        a.best_front_elem = frontel
    elseif hat < 0 # candidate replaced frontier element
        # update best_front_elem if it was replaced or if the new aggscore is better
        # NOTE: when the old best element is replaced, the new aggscore might be
        #       worse due to differences in aggscore and eps-progress metrics
        d = a.best_front_elem.fitness.agg - frontel.fitness.agg
        if  a.best_front_elem.fitness.index == frontel.fitness.index ||
           (d > zero(d) && is_minimizing(a.fit_scheme)) ||
           (d < zero(d) && !is_minimizing(a.fit_scheme))
            #@debug "New best candidate"
            a.best_front_elem = frontel
        end
    end
    if length(a.frontier) <= a.max_size
        a.n_oversize_inserts = 0 # reset the counter since the size is ok
    end
    return a
end

# actually this methods should never be called because the fitness
# is already indexes within the method
add_candidate!(a::EpsBoxArchive{N,F}, cand_fitness::NTuple{N,F},
               candidate::AbstractIndividual, tag::Int=0, num_fevals::Int=-1) where {N,F} =
    add_candidate!(a, archived_fitness(cand_fitness, a), candidate, tag, num_fevals)

# called by check_stop_condition(e::Evaluator, ctrl)
function check_stop_condition(a::EpsBoxArchive, p::OptimizationProblem, ctrl)
    if ctrl.max_steps_without_progress > 0 &&
        noprogress_streak(a, since_restart=false) > ctrl.max_steps_without_progress
        return "No epsilon-progress for more than $(ctrl.max_steps_without_progress) iterations"
    elseif a.n_oversize_inserts >= 10
        # notify that the last 10 inserts were to the oversized archive
        # that means that the ϵ quantization steps are too small
        return "Pareto frontier size ($(length(a.frontier))) exceeded maximum ($(a.max_size))"
    else
        return ""
    end
end
