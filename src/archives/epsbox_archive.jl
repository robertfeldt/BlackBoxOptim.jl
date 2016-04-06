"""
    Individual stored in `EpsBoxArchive`.
"""
immutable EpsBoxFrontierIndividual{N,F<:Number} <: ArchivedIndividual{IndexedTupleFitness{N,F}}
    fitness::IndexedTupleFitness{N,F}
    params::Individual
    tag::Int                            # tag of the individual (e.g. gen.op. ID)
    timestamp::Float64                  # when archived
    num_fevals::Int                     # number of fitness evaluations so far

    Base.call{N,F}(::Type{EpsBoxFrontierIndividual}, fitness::IndexedTupleFitness{N,F},
                   params, tag, num_fevals) =
        new{N,F}(fitness, deepcopy(params), tag, time(), num_fevals)
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
  candidates_without_progress::Int  # last Number of calls to add_candidate!() without ϵ-progress

  max_size::Int         # maximal frontier size
  # TODO allow different frontier containers?
  # see e.g. Altwaijry & Menai "Data Structures in Multi-Objective Evolutionary Algorithms", 2012
  frontier::Vector{EpsBoxFrontierIndividual{N,F}}  # candidates along the fitness Pareto frontier

  function Base.call{N,F}(::Type{EpsBoxArchive}, fit_scheme::EpsBoxDominanceFitnessScheme{N,F}; max_size::Integer = 1_000_000)
    new{N,F,typeof(fit_scheme)}(fit_scheme, time(), 0, 0, 0, max_size, EpsBoxFrontierIndividual{N,F}[])
  end

  Base.call{N,F}(::Type{EpsBoxArchive}, fit_scheme::EpsBoxDominanceFitnessScheme{N,F}, params::Parameters) =
    EpsBoxArchive(fit_scheme, max_size=params[:MaxArchiveSize])
end

const EpsBoxArchive_DefaultParameters = ParamsDict(
  :MaxArchiveSize => 10_000,
)

fitness_scheme(a::EpsBoxArchive) = a.fit_scheme
# EpsBoxArchive stores indexed fitness
Base.length(a::EpsBoxArchive) = length(a.frontier)
Base.isempty(a::EpsBoxArchive) = isempty(a.frontier)
capacity(a::EpsBoxArchive) = a.max_size
numdims(a::EpsBoxArchive) = !isempty(a.frontier) ? length(a.frontier[1].params) : 0
pareto_frontier(a::EpsBoxArchive) = a.frontier

"""
    `candidates_without_progress(a::EpsBoxArchive)`

    Get the number of `add_candidate!()` calls since the last ϵ-progress.
"""
candidates_without_progress(a::EpsBoxArchive) = a.candidates_without_progress

Base.getindex(a::EpsBoxArchive, i::Integer) = a.frontier[i]

best_candidate(a::EpsBoxArchive) = a.frontier[a.best_candidate_ix].params
best_fitness(a::EpsBoxArchive) = a.best_candidate_ix > 0 ? fitness(a.frontier[a.best_candidate_ix]) : nafitness(fitness_scheme(a))

"""
    `tagcounts(a::EpsBoxArchive)`

    Count the tags of individuals on the ϵ-box frontier.
    Returns the `tag`→`count` dictionary.
"""
function tagcounts(a::EpsBoxArchive)
    res = Dict{Int,Int}()
    for i in eachindex(a.frontier)
        curtag = tag(a.frontier[i])
        curcounts = get!(res, curtag, 0)
        res[curtag] = curcounts+1
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

    has_progress = true
    updated_frontier_ix = 0
    for i in length(a.frontier):-1:1
        hat, index_match = hat_compare(cand_fitness, fitness(a.frontier[i]), a.fit_scheme)
        has_progress = has_progress && !index_match
        if hat < 0
            # new fitness dominates one in the archive
            if updated_frontier_ix==0
                #info("Replaced the dominated element $i on the frontier")
                # replace the fitness to minimize memory operations
                # note - if i was the best candidate index, it stays
                a.frontier[i] = EpsBoxFrontierIndividual(cand_fitness, candidate, tag, num_fevals)
                updated_frontier_ix = i
            else
                # already replaced the other dominated fitness, delete this one
                #info("Removed the dominated element $i from the frontier")
                splice!(a.frontier, i)
                if updated_frontier_ix > i
                    updated_frontier_ix -= 1
                end
                if a.best_candidate_ix > i
                    a.best_candidate_ix -= 1
                elseif a.best_candidate_ix == i
                    #info("New best candidate")
                    a.best_candidate_ix = updated_frontier_ix
                end
            end
        elseif hat > 0 || (hat == 0 && index_match &&
                           isapprox(cand_fitness.agg, fitness(a.frontier[i]).agg, rtol=0.0, atol=100eps(Float64)) &&
                           isapprox(cand_fitness.dist, fitness(a.frontier[i]).dist, rtol=0.0, atol=100eps(Float64)))
            # the new fitness dominates or is equal to the element in the archive, ignore
            has_progress = false
            updated_frontier_ix = -1 # just something nonzero to prevent appending
            break
        end
    end
    if updated_frontier_ix == 0 # non-dominated candidate, append to the frontier
        #info("Appended non-dominated element to the frontier")
        push!(a.frontier, EpsBoxFrontierIndividual(cand_fitness, candidate, tag, num_fevals))
        if length(a.frontier) > a.max_size
            throw(error("Pareto frontier exceeds maximum size"))
        end
        updated_frontier_ix = length(a.frontier)
        # check if the new candidate has better aggregate score
        if a.best_candidate_ix==0
            a.best_candidate_ix = updated_frontier_ix
        else
            d = a.frontier[a.best_candidate_ix].fitness.agg - a.frontier[updated_frontier_ix].fitness.agg
            if (d > zero(d) && is_minimizing(a.fit_scheme)) || (d < zero(d) && !is_minimizing(a.fit_scheme))
                #info("New best candidate")
                a.best_candidate_ix = updated_frontier_ix
            end
        end
    end
    if has_progress # non-dominated solution that has some indices different from the existing ones
        a.candidates_without_progress = 0
    else
        a.candidates_without_progress += 1
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
check_stop_condition(a::EpsBoxArchive, p::OptimizationProblem, ctrl) = ""
