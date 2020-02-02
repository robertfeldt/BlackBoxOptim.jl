"""
Internal data for the worker process of the parallel evaluator.
"""
struct ParallelEvaluatorWorker{P<:OptimizationProblem}
    problem::P

    ParallelEvaluatorWorker(problem::P) where {P<:OptimizationProblem} =
        new{P}(problem)
end

fitness(params::Individual, worker::ParallelEvaluatorWorker) =
    fitness(params, worker.problem)

const ChannelRef{T} = RemoteChannel{Channel{T}}
const ParallelEvaluatorWorkerRef{P} = ChannelRef{ParallelEvaluatorWorker{P}}

fitness(params::Individual, worker_ref::ParallelEvaluatorWorkerRef{P}) where P =
    fitness(params, fetch(fetch(worker_ref))::ParallelEvaluatorWorker{P})

"""
Current state of fitness function evaluation for the vector of candidates.
"""
mutable struct ParallelEvaluationState{FA, A}
    archive::A
    candidates::Vector{Candidate{FA}}  # candidates to calculate fitness for

    worker_busy::Vector{Bool}         # if the worker is busy calculating
    retry_queue::Vector{Int}          # queue to candidate indices to retry

    worker_finished::Condition        # gets notified each time worker is done
                                      # fitness calculation
    next_index::Int                   # index of the next candidate to calculate fitness

    function ParallelEvaluationState(archive::A, nworkers::Int) where {A<:Archive}
        FA = fitness_type(archive)
        new{FA,A}(archive, Vector{Candidate{FA}}(),
                  fill(false, nworkers), Vector{Int}(), Condition(), 0)
    end
end

is_stopped(estate::ParallelEvaluationState) = estate.next_index == 0

abort!(estate::ParallelEvaluationState) = (estate.next_index = 0)

"""
Reset the current `ParallelEvaluationState` and the vector
of candidates that need fitness evaluation.
"""
function reset!(estate::ParallelEvaluationState{FA}, candidates::Vector{Candidate{FA}}) where FA
    estate.candidates = candidates
    fill!(estate.worker_busy, false)
    empty!(estate.retry_queue)
    estate.next_index = 1
    return estate
end

"""
Get the index of the next candidate for evaluation
based on the Base.pmap() code.
"""
function next_candidate!(estate::ParallelEvaluationState, worker_ix::Int)
    @assert !estate.worker_busy[worker_ix]

    task_ix = 0
    if !is_stopped(estate) && estate.next_index <= length(estate.candidates)
        # find the next candidate with unevaluated fitness
        while !is_stopped(estate) && estate.next_index <= length(estate.candidates)
            if isnafitness(estate.candidates[estate.next_index].fitness, fitness_scheme(estate.archive))
                task_ix = estate.next_index
                estate.next_index += 1
                break;
            else
                estate.next_index += 1
            end
        end
    end
    if task_ix == 0
        if !isempty(estate.retry_queue)
            task_ix = popfirst!(estate.retry_queue)
        else
            # Handles the condition where we have finished processing the requested
            # lsts as well as any retryqueue entries, but there are still some jobs
            # active that may result in an error and have to be retried.
            while any(estate.worker_busy)
                wait(estate.worker_finished)
                if !isempty(estate.retry_queue)
                    task_ix = popfirst!(estate.retry_queue)
                    break
                end
            end
        end
    end
    estate.worker_busy[worker_ix] = task_ix != 0
    return task_ix
end

"""
Notify that the worker process is finished and reset its busy flag.
"""
function worker_finished!(estate::ParallelEvaluationState, worker_ix::Int)
    @assert estate.worker_busy[worker_ix]
    estate.worker_busy[worker_ix] = false
    notify(estate.worker_finished; all=true)
    return estate
end

"""
Fitness evaluator that distributes candidates fitness calculation
among several worker processes.
"""
mutable struct ParallelEvaluator{F, FA, FS, P<:OptimizationProblem, A<:Archive} <: Evaluator{P}
    problem::P
    archive::A
    num_evals::Int
    last_fitness::F
    arch_nafitness::FA  # NA fitness

    worker_refs::Vector{ParallelEvaluatorWorkerRef{P}}
    eval_state::ParallelEvaluationState{FA, A}
end

function ParallelEvaluator(
        problem::P, archive::A;
        pids::Vector{Int} = workers()) where {P<:OptimizationProblem, A<:Archive}
    fs = fitness_scheme(problem)
    F = fitness_type(fs)
    FA = fitness_type(archive)
    ParallelEvaluator{F, FA, typeof(fs), P, A}(problem,
        archive,
        0, nafitness(fs), nafitness(FA),
        [RemoteChannel(function ()
                     # create fake channel and put problem there
                     ch = Channel{ParallelEvaluatorWorker{P}}(1)
                     put!(ch, ParallelEvaluatorWorker(copy(problem)))
                     ch
                   end, pid) for pid in pids],
        ParallelEvaluationState(archive, length(pids))
    )
end

ParallelEvaluator(problem::OptimizationProblem;
        pids::Vector{Int} = workers(),
        archiveCapacity::Integer = 10) =
    ParallelEvaluator(problem, TopListArchive(fitness_scheme(problem), numdims(problem), archiveCapacity),
                      pids=pids)

num_evals(e::ParallelEvaluator) = e.num_evals

function update_fitness!(f::Any, e::ParallelEvaluator, candidates::Any;
                         force::Bool=false)
    # FIXME use force
    reset!(e.eval_state, candidates)

    # based on pmap() code
    @sync begin
        for (widx, wref) in enumerate(e.worker_refs)
            @async begin
                while (candi_ix = next_candidate!(e.eval_state, widx)) != 0
                    try
                        candi = candidates[candi_ix]
                        candi_fitness = remotecall_fetch(fitness, wref.where, candi.params, wref)
                        candi.fitness = archived_fitness(candi_fitness, e.archive)
                        e.last_fitness = candi_fitness
                        e.num_evals += 1
                        add_candidate!(e.archive, candi.fitness, candi.params, e.num_evals)
                        (f !== nothing) && f(candi)
                    catch ex
                        abort!(e.eval_state) # when one worker fails, the whole process is aborted
                        rethrow(ex)
                    finally
                        worker_finished!(e.eval_state, widx)
                    end
                end
            end
        end
    end

    return candidates
end

# FIXME it's not efficient to calculate fitness like that with `ParallelEvaluator`
update_fitness!(f::Any, e::ParallelEvaluator, candidate::Candidate;
                force::Bool=false) = update_fitness!(f, e, [candidate])

# FIXME it's not efficient to calculate fitness like that with `ParallelEvaluator`
function fitness(params::Individual, e::ParallelEvaluator{F, FA}, tag::Int=0) where {F, FA}
    candi = Candidate{FA}(params, -1, e.arch_nafitness, nothing, tag)
    update_fitness!(e, [candi])
    return candi.fitness
end
