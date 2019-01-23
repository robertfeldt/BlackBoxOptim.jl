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
mutable struct ParallelEvaluationState{F, FS}
    fitness_scheme::FS
    candidates::Vector{Candidate{F}}  # candidates to calculate fitness for

    worker_busy::Vector{Bool}         # if the worker is busy calculating
    retry_queue::Vector{Int}          # queue to candidate indices to retry

    worker_finished::Condition        # gets notified each time worker is done
                                      # fitness calculation
    next_index::Int                   # index of the next candidate to calculate fitness

    function ParallelEvaluationState(fitness_scheme::FS,
                                     nworkers::Int) where {FS<:FitnessScheme}
        F = fitness_type(fitness_scheme)
        new{F,FS}(fitness_scheme, Vector{Candidate{F}}(),
                  fill(false, nworkers), Vector{Int}(), Condition(), 0)
    end
end

is_stopped(estate::ParallelEvaluationState) = estate.next_index == 0

abort!(estate::ParallelEvaluationState) = (estate.next_index = 0)

"""
Reset the current `ParallelEvaluationState` and the vector
of candidates that need fitness evaluation.
"""
function reset!(estate::ParallelEvaluationState{F}, candidates::Vector{Candidate{F}}) where F
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
            if isnafitness(estate.candidates[estate.next_index].fitness, estate.fitness_scheme)
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
mutable struct ParallelEvaluator{F, FS, P<:OptimizationProblem} <: Evaluator{P}
    problem::P
    archive::Archive
    num_evals::Int
    last_fitness::F

    worker_refs::Vector{ParallelEvaluatorWorkerRef{P}}
    eval_state::ParallelEvaluationState{F, FS}
end

function ParallelEvaluator(
        problem::P, archive::Archive;
        pids::Vector{Int} = workers()) where {P<:OptimizationProblem}
    fs = fitness_scheme(problem)
    ParallelEvaluator{fitness_type(fs), typeof(fs), P}(problem,
        archive,
        0, nafitness(fs),
        [RemoteChannel(function ()
                     # create fake channel and put problem there
                     ch = Channel{ParallelEvaluatorWorker{P}}(1)
                     put!(ch, ParallelEvaluatorWorker(copy(problem)))
                     ch
                   end, pid) for pid in pids],
        ParallelEvaluationState(fs, length(pids))
    )
end

ParallelEvaluator(problem::OptimizationProblem;
        pids::Vector{Int} = workers(),
        archiveCapacity::Integer = 10) =
    ParallelEvaluator(problem, TopListArchive(fitness_scheme(problem), numdims(problem), archiveCapacity),
                      pids=pids)

num_evals(e::ParallelEvaluator) = e.num_evals

function update_fitness!(e::ParallelEvaluator{F}, candidates::Vector{Candidate{F}}) where F
    reset!(e.eval_state, candidates)

    # based on pmap() code
    @sync begin
        for (widx, wref) in enumerate(e.worker_refs)
            @async begin
                while (candi_ix = next_candidate!(e.eval_state, widx)) != 0
                    try
                        candi = candidates[candi_ix]
                        candi.fitness = remotecall_fetch(fitness, wref.where, candi.params, wref)
                        e.last_fitness = candi.fitness
                        e.num_evals += 1
                        add_candidate!(e.archive, candi.fitness, candi.params, e.num_evals)
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
function fitness(params::Individual, e::ParallelEvaluator{F}) where F
    candi = Candidate{F}(params)
    update_fitness!(e, [candi])
    return candi.fitness
end
