"""
    Individual stored in `EpsBoxArchive`.
"""
immutable EpsBoxFrontierIndividual{N,F<:Number} <: ArchivedIndividual{IndexedTupleFitness{N,F}}
    fitness::IndexedTupleFitness{N,F}
    params::Individual
    tag::Int                            # tag of the individual (e.g. gen.op. ID)
    timestamp::Float64                  # when archived
    num_fevals::Int                     # number of fitness evaluations so far
    n_restarts::Int                     # the number of method restarts so far

    @compat (::Type{EpsBoxFrontierIndividual}){F,N}(fitness::IndexedTupleFitness{N,F},
                   params, tag, num_fevals, n_restarts) =
        new{N,F}(fitness, params, tag, time(), num_fevals, n_restarts)
end

"""
    ϵ-box archive saves only the solutions that are not ϵ-box
    dominated by any other solutions in the archive.

    It also counts the number of candidate solutions that have been added
    and how many ϵ-box progresses have been made.
"""
type EpsBoxArchive{N,F,FS<:EpsBoxDominanceFitnessScheme} <: Archive{NTuple{N,F},IndexedTupleFitness{N,F},FS}
  fit_scheme::FS        # Fitness scheme used
  start_time::Float64   # Time when archive created, we use this to approximate the starting time for the opt...

  num_candidates::Int               # Number of calls to add_candidate!()
  best_candidate_ix::Int            # the index of the candidate with the best aggregated fitness
  last_progress::Int                # when (wrt num_candidates) last ϵ-progress has occured
  last_restart::Int                 # when (wrt num_dlast) last restart has occured
  n_restarts::Int                   # the counter of the method restarts

  len::Int              # current frontier size
  max_size::Int         # maximal frontier size
  # TODO allow different frontier containers?
  # see e.g. Altwaijry & Menai "Data Structures in Multi-Objective Evolutionary Algorithms", 2012
  frontier::Vector{EpsBoxFrontierIndividual{N,F}}  # candidates along the fitness Pareto frontier
  frontier_isoccupied::BitVector # true if given frontier element is occupied

  @compat function (::Type{EpsBoxArchive}){N,F}(fit_scheme::EpsBoxDominanceFitnessScheme{N,F}; max_size::Integer = 1_000_000)
    new{N,F,typeof(fit_scheme)}(fit_scheme, time(), 0, 0, 0, 0, 0, 0, max_size,
                                sizehint!(Vector{EpsBoxFrontierIndividual{N,F}}(), 64),
                                sizehint!(BitVector(), 64))
  end

  @compat (::Type{EpsBoxArchive}, fit_scheme::EpsBoxDominanceFitnessScheme{N,F}, params::Parameters) =
    EpsBoxArchive(fit_scheme, max_size=params[:MaxArchiveSize])
end

const EpsBoxArchive_DefaultParameters = ParamsDict(
  :MaxArchiveSize => 10_000,
)

fitness_scheme(a::EpsBoxArchive) = a.fit_scheme
# EpsBoxArchive stores indexed fitness
Base.length(a::EpsBoxArchive) = a.len
Base.isempty(a::EpsBoxArchive) = a.len == 0
capacity(a::EpsBoxArchive) = a.max_size
numdims(a::EpsBoxArchive) = !isempty(a.frontier) ? length(a.frontier[1].params) : 0

"""
    Iterates occupied elements of the `archive.frontier`.
"""
immutable EpsBoxArchiveFrontierIterator{A<:EpsBoxArchive}
    archive::A
    @compat (::Type{EpsBoxArchiveFrontierIterator}){A<:EpsBoxArchive}(a::A) = new{A}(a)
end

@inline Base.length(it::EpsBoxArchiveFrontierIterator) = length(it.archive)
@inline Base.start(it::EpsBoxArchiveFrontierIterator) = findfirst(it.archive.frontier_isoccupied)
@inline Base.done(it::EpsBoxArchiveFrontierIterator, ix::Integer) = ix == 0
@inline Base.next(it::EpsBoxArchiveFrontierIterator, ix::Integer) = (it.archive.frontier[ix], findnext(it.archive.frontier_isoccupied, ix+1))

"""
    Get the iterator to the individuals on the Pareto frontier.
"""
pareto_frontier(a::EpsBoxArchive) = EpsBoxArchiveFrontierIterator(a)

function occupied_frontier_indices(a::EpsBoxArchive)
    ixs = sizehint!(Vector{Int}(), a.len)
    i = findfirst(a.frontier_isoccupied)
    while i > 0
        push!(ixs, i)
        i = findnext(a.frontier_isoccupied, i+1)
    end
    @assert length(ixs) == a.len
    return ixs
end

"""
    Get random occupied Pareto frontier index.
    Returns 0 if frontier is empty.
"""
function rand_frontier_index(a::EpsBoxArchive)
    if a.len == 0
        return 0
    end
    # narrow the range of random indices
    il = findfirst(a.frontier_isoccupied)
    iu = findlast(a.frontier_isoccupied)
    # generate random indices until occupied one found
    # FIXME any better scheme to get uniform samples from occupied indices?
    i = rand(il:iu)
    while !a.frontier_isoccupied[i]
        i = rand(il:iu)
    end
    return i
end

"""
    `noprogress_streak(a::EpsBoxArchive, [since_restart])`

    Get the number of `add_candidate!()` calls since the last ϵ-progress.
    If `since_restart` is specified, the number is relative to the last
    restart.
"""
noprogress_streak(a::EpsBoxArchive; since_restart::Bool=false) =
    since_restart ?
        a.num_candidates - max(a.last_progress, a.last_restart) :
        a.num_candidates - a.last_progress

best_candidate(a::EpsBoxArchive) = a.frontier[a.best_candidate_ix].params
best_fitness(a::EpsBoxArchive) = a.best_candidate_ix > 0 ? fitness(a.frontier[a.best_candidate_ix]) : nafitness(fitness_scheme(a))
function notify!(a::EpsBoxArchive, event::Symbol)
    if event == :restart
        a.n_restarts += 1
        a.last_restart = a.num_candidates
        # TODO check if the pareto frontier needs to be compactified
    end
    a
end

"""
    `tagcounts(a::EpsBoxArchive)`

    Count the tags of individuals on the ϵ-box frontier.
    Each restart the individual remains in the frontier discounts it by `θ`.

    Returns the `tag`→`count` dictionary.
"""
function tagcounts(a::EpsBoxArchive, θ::Number = 1.0)
    (0.0 < θ <= 1.0) || throw(ArgumentError("θ ($θ) should be in (0.0, 1.0] range"))
    res = Dict{Int,Float64}()
    i = findfirst(a.frontier_isoccupied)
    while i > 0
        curtag = tag(a.frontier[i])
        if curtag > 0
            curcounts = get!(res, curtag, 0.0)
            res[curtag] = curcounts+θ^(a.n_restarts-a.frontier[i].n_restarts)
        end
        i = findnext(a.frontier_isoccupied, i+1)
    end
    return res
end

function add_candidate!{N,F}(a::EpsBoxArchive{N,F}, cand_fitness::IndexedTupleFitness{N,F},
                             candidate, tag::Int=0, num_fevals::Int=-1)
    a.num_candidates += 1
    if num_fevals == -1
        num_fevals = a.num_candidates
    end
    #info("New fitness: ", cand_fitness.orig, " agg=", cand_fitness.agg)
    #info("Params: ", candidate)

    has_progress = false
    updated_frontier_ix = 0
    i = findfirst(a.frontier_isoccupied)
    while i > 0
        front_fitness = fitness(a.frontier[i])
        hat, index_match = hat_compare(cand_fitness, front_fitness, a.fit_scheme)
        if hat < 0
            # new fitness dominates the one in the archive
            if updated_frontier_ix==0
                #info("Replaced the dominated element $i on the frontier")
                # replace the fitness to minimize memory operations
                # note - if i was the best candidate index, it stays
                new_params = copy!(a.frontier[i].params, candidate)
                a.frontier[i] = EpsBoxFrontierIndividual(cand_fitness, new_params, tag, num_fevals, a.n_restarts)
                updated_frontier_ix = i
                if index_match # we replace the element with the same index, so the domination pattern to the other elements should not change
                    break
                else
                    has_progress = true
                end
            else
                # already replaced the other dominated fitness, exclude this one from the frontier
                #info("Removed the dominated element $i from the frontier")
                @assert !index_match # if matching, it means that the old frontier[updated_frontier_ix] _was_ dominated by frontier[i]
                a.frontier_isoccupied[i] = false
                a.len -= 1
                if a.best_candidate_ix == i
                    #info("New best candidate")
                    a.best_candidate_ix = updated_frontier_ix
                end
            end
        elseif hat > 0 || (hat == 0 && index_match)
            # the candidate fitness is worse or just the same as in the frontier, don't insert it
            @assert updated_frontier_ix == 0 # should not have been inserted into the frontier before
            updated_frontier_ix = -1 # just something nonzero to prevent appending
            break
        end
        i = findnext(a.frontier_isoccupied, i+1)
    end
    if updated_frontier_ix == 0 # non-dominated candidate, append to the frontier
        updated_frontier_ix = findfirst(a.frontier_isoccupied, false) # first unoccupied
        if updated_frontier_ix > 0
            # replace the deactivated frontier element and activate it
            new_params = copy!(a.frontier[updated_frontier_ix].params, candidate)
            a.frontier[updated_frontier_ix] = EpsBoxFrontierIndividual(cand_fitness, new_params, tag, num_fevals, a.n_restarts)
            a.frontier_isoccupied[updated_frontier_ix] = true
        else
            # all frontier elements are occupied, append another one
            push!(a.frontier, EpsBoxFrontierIndividual(cand_fitness, deepcopy(candidate), tag, num_fevals, a.n_restarts))
            push!(a.frontier_isoccupied, true)
            updated_frontier_ix = length(a.frontier)
        end
        a.len += 1
        has_progress = true # no existing element with the same index, otherwise it would have been replaced
        #info("Appended non-dominated element to the frontier")
        if length(a.frontier) > a.max_size
            throw(error("Pareto frontier exceeds maximum size"))
        end
    end
    if updated_frontier_ix > 0
        # check if the new candidate has better aggregate score
        if a.best_candidate_ix==0
            a.best_candidate_ix = updated_frontier_ix
        elseif a.best_candidate_ix != updated_frontier_ix
            d = a.frontier[a.best_candidate_ix].fitness.agg - cand_fitness.agg
            if (d > zero(d) && is_minimizing(a.fit_scheme)) || (d < zero(d) && !is_minimizing(a.fit_scheme))
                #info("New best candidate")
                a.best_candidate_ix = updated_frontier_ix
            end
        end
    end
    if has_progress # non-dominated solution that has some eps-indices different from the existing ones
        a.last_progress = a.num_candidates
    end
    return a
end

archived_fitness{N,F}(fitness::NTuple{N,F}, a::EpsBoxArchive{N,F}) =
    convert(IndexedTupleFitness, fitness, a.fit_scheme)

# actually this methods should never be called because the fitness
# is already indexes within the method
add_candidate!{N,F}(a::EpsBoxArchive{N,F}, cand_fitness::NTuple{N,F},
                    candidate, tag::Int=0, num_fevals::Int=-1) =
    add_candidate!(a, archived_fitness(cand_fitness, a), candidate, tag, num_fevals)

# called by check_stop_condition(e::Evaluator, ctrl)
function check_stop_condition(a::EpsBoxArchive, p::OptimizationProblem, ctrl)
    ctrl.max_steps_without_progress > 0 && noprogress_streak(a, since_restart=false) > ctrl.max_steps_without_progress ?
        "No epsilon-progress for more than $(ctrl.max_steps_without_progress) iterations" : ""
end
