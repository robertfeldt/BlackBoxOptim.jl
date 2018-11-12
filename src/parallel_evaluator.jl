using Distributed, SharedArrays

#= master <-> worker params_status/fitness_status codes =#

const PEStatus_OK = 0 # message received; positive statuses = job_id
const PEStatus_Shutdown = -1
const PEStatus_Error = -2

"""
Internal data for the worker process of the parallel evaluator.
"""
mutable struct ParallelEvaluatorWorker{T, P<:OptimizationProblem}
    id::Int                             # worker ID
    problem::P
    param_status::SharedVector{Int}     # master notifies worker about new requests
    shared_param::SharedVector{T}       # master puts candidates parameters
    fitness_status::SharedVector{Int}   # worker notifies master about completed evaluation
    shared_fitness::SharedVector{T}     # worker puts calculated fitness

    ParallelEvaluatorWorker(
        id::Int, problem::P,
        param_status::SharedVector{Int}, shared_param::SharedVector{T},
        fitness_status::SharedVector{Int}, shared_fitness::SharedVector{T}
    ) where {T, P<:OptimizationProblem} =
        new{T,P}(id, problem, param_status, shared_param, fitness_status, shared_fitness)
end

param_status(worker::ParallelEvaluatorWorker) = @inbounds first(worker.param_status)
fitness_status(worker::ParallelEvaluatorWorker) = @inbounds first(worker.fitness_status)

# run the wrapper (called in the "main" task)
function run!(worker::ParallelEvaluatorWorker)
    while true
        #@info "Checking in worker #$(worker.id)"
        # continuously poll the worker for the delivery notification for
        # the last job or for the new job notification
        i = 0
        while param_status(worker) == PEStatus_OK || fitness_status(worker) > 0
            if (i+=1) > 1000 # allow executing the other tasks once in a while
                yield()
                i = 0
            end
        end
        # process the new worker status
        p_status = param_status(worker)
        if p_status == PEStatus_Shutdown # master stopping
            worker.fitness_status[1] = PEStatus_Shutdown # notify worker stopped
            break # shutdown!() called
        elseif p_status == PEStatus_Error # master error (currently not set?)
            worker.fitness_status[1] = PEStatus_Shutdown # stopped after receiving an error
            break
        elseif p_status > 0 # new job (status is job_id)
            worker.param_status[1] = PEStatus_OK # received, reset the statuts
            #@info "PE worker #$(worker.id): got job"
            @inbounds setfitness!(worker.shared_fitness, fitness(worker.shared_param, worker.problem))
            worker.fitness_status[1] = p_status # fitness ready
        end
    end
end

"""
Create and run the evaluator worker.
The function that the master process spawns at each worker process.
"""
function run_parallel_evaluator_worker(id::Int,
                    worker_ready::RemoteChannel{Channel{Int}},
                    problem::OptimizationProblem,
                    param_status::SharedVector{Int},
                    shared_param::SharedVector{Float64},
                    fitness_status::SharedVector{Int},
                    shared_fitness::SharedVector{Float64})
    @info "Initializing ParallelEvaluator worker #$id at task=$(myid())"
    worker = nothing
    try
        worker = ParallelEvaluatorWorker(id, deepcopy(problem),
                     param_status, shared_param,
                     fitness_status, shared_fitness)
    catch e
        # send -id to notify about an error and to release
        # the master from waiting for worker readiness
        @warn "Exception at ParallelEvaluatorWorker initialization: $e"
        put!(worker_ready, -id)
        rethrow(e)
    end
    # create immigrants receiving tasks=#
    put!(worker_ready, id)
    @info "Running worker #$id..."
    try
        run!(worker)
    catch e
        # send error candidate to notify about an error and to release
        # the master from waiting for worker messages
        @warn "Exception while running ParallelEvaluatorWorker: $e"
        worker.fitness_status[1] = PEStatus_Error
        rethrow(e)
    end
    @info "Worker #$id stopped"
    nothing
end

const PECandidateDict{FA} = Dict{Int, Candidate{FA}}

"""
Fitness evaluator that asynchronously distributes calculation
among several worker processes.

Currently the overhead of coordinating parallel processes is relatively high,
so it's recommended to use ParallelEvaluator only when fitness calculation takes
considerable time.
"""
mutable struct ParallelEvaluator{F, FA, T, FS, P<:OptimizationProblem, A<:Archive} <: Evaluator{P}
    problem::P          # optimization problem
    archive::A          # archive where good candidates are automatically stored
    num_evals::Int      # fitness evaluations counter
    last_fitness::F     # last fitness
    arch_nafitness::FA  # NA fitness

    params_status::Vector{SharedVector{Int}}    # master notifies workers about new requests
    shared_params::Vector{SharedVector{T}}      # master puts candidates parameters
    fitnesses_status::Vector{SharedVector{Int}} # workers notify master about completed evaluation
    shared_fitnesses::Vector{SharedVector{T}}   # workers put calculated fitness

    fitness_slots::Base.Semaphore               # gets acquired when a worker needs to be assigned to a task;
                                                # used to organize waiting when all workers are busy

    waiting_candidates::PECandidateDict{FA}     # candidates waiting their fitness calculation to be completed
    unclaimed_candidates::PECandidateDict{FA}   # candidates with calculated fitness that were not yet checked for completion (by process_completed())

    busy_workers::BitVector         # workers available for accepting jobs
    job_assignment::ReentrantLock   # lock to provide exclusive access to busy_workers

    max_seq_done_job::Int   # all jobs from 1 to max_seq_done_job are done
    max_done_job::Int       # max Id of done job
    done_jobs::BitSet       # done job Ids beyond max_seq_done_job

    next_job_id::Int        # ID to assign for the next job

    is_stopping::Bool       # whether the evaluator is in shutdown sequence

    worker_refs::Vector{Future} # references to worker processes
    workers_listener::Task  # task in the main process that runs workers_listener!()

    function ParallelEvaluator(
        problem::P, archive::A;
        pids::AbstractVector{Int} = workers()
    ) where {P<:OptimizationProblem, A<:Archive}
        fs = fitness_scheme(problem)
        F = fitness_type(fs)
        T = fitness_eltype(fs)
        FA = fitness_type(archive)

        etor = new{F, FA, T, typeof(fs), P, A}(
            problem, archive,
            0, nafitness(fs), nafitness(FA),
            [fill!(SharedArray{Int}((2,), pids=vcat(pid,[myid()])), 0) for pid in pids],
            [SharedArray{T}((numdims(problem),), pids=vcat(pid,[myid()])) for pid in pids],
            [fill!(SharedArray{Int}((2,), pids=vcat(pid,[myid()])), 0) for pid in pids],
            [SharedArray{T}((numobjectives(fs),), pids=vcat(pid,[myid()])) for pid in pids],
            Base.Semaphore(length(pids)),
            PECandidateDict{FA}(), PECandidateDict{FA}(),
            falses(length(pids)), ReentrantLock(), 0, 0, BitSet(),
            1, false
        )
        etor.worker_refs = _create_workers(etor, pids)
        etor.workers_listener = @async workers_listener!(etor)

        #finalizer(etor, _shutdown!)
        return etor
    end
end

ParallelEvaluator(
    problem::OptimizationProblem;
    pids::AbstractVector{Int} = workers(),
    archiveCapacity::Integer = 10) =
    ParallelEvaluator(problem, TopListArchive(fitness_scheme(problem), numdims(problem), archiveCapacity),
                      pids=pids)

nworkers(etor::ParallelEvaluator) = length(etor.worker_refs)
queue_capacity(etor::ParallelEvaluator) = nworkers(etor)

"""
Count the candidates submitted (including the completed ones),
but not yet claimed.
"""
queue_length(etor::ParallelEvaluator) = length(etor.waiting_candidates) + length(etor.unclaimed_candidates)

num_evals(etor::ParallelEvaluator) = etor.num_evals

is_stopping(etor::ParallelEvaluator) = etor.is_stopping

# check that worker is stil running.
# If running, its RemoteChannels should not be ready,
# but if there was exception in the worker,
# it would be thrown into the main thread
function check_worker_running(worker::Future)
    if isready(worker)
        worker_res = fetch(worker) # fetch the worker, this should trigger an exception
        # no exception, but the worker should not be ready
        error("Worker at pid=$(worker.where) has finished before the master shutdown: $worker_res")
    end
    return true
end

function _create_workers(etor::ParallelEvaluator, pids::AbstractVector{Int})
    @info "Initializing parallel workers..."
    workers_ready = RemoteChannel(() -> Channel{Int}(length(pids))) # FIXME do we need to wait for the worker?

    # spawn workers
    problem = etor.problem
    params_status = etor.params_status
    shared_params = etor.shared_params
    fitnesses_status = etor.fitnesses_status
    shared_fitnesses = etor.shared_fitnesses

    worker_refs = Future[
        @spawnat(pid, run_parallel_evaluator_worker(i, workers_ready, problem,
                        params_status[i], shared_params[i],
                        fitnesses_status[i], shared_fitnesses[i]))
        for (i, pid) in enumerate(pids)]
    #@assert !isready(ppopt.is_started)
    # wait until all the workers are started
    @info "Waiting for the workers to be ready..."
    # FIXME is it required?
    nready = 0
    while nready < length(pids)
        check_worker_running.(worker_refs)
        worker_id = take!(workers_ready)
        if worker_id < 0
            # worker failed to initialize, reading its task would throw an exception
            check_worker_running(worker_refs[-worker_id])
            error("Exception in the worker #$(-worker_id), but all workers still running")
        end
        @info "  Worker #$worker_id is ready"
        nready += 1
    end
    @info "All workers ready"
    return worker_refs
end

function shutdown!(etor::ParallelEvaluator)
    @info "shutdown!(ParallelEvaluator)"
    etor.is_stopping && error("Cannot shutdown!(ParallelEvaluator) twice")
    etor.is_stopping = true
    # notify the workers that they should shutdown (each worker should pick exactly one message)
    _shutdown!(etor)
    # resume workers listener if it is waiting for the new jobs
    lock(etor.job_assignment)
    unlock(etor.job_assignment)
    # wait for all the workers
    for i in 1:nworkers(etor)
        Base.acquire(etor.fitness_slots)
    end
    @assert !any(etor.busy_workers) "Some workers have not finished"
    # release any waiting
    for i in 1:nworkers(etor)
        Base.release(etor.fitness_slots)
    end
end

function _shutdown!(etor::ParallelEvaluator)
    #@info "_shutdown!(ParallelEvaluator)"
    if !etor.is_stopping
        etor.is_stopping = true
        #close(etor.in_fitnesses)
        #close(etor.out_individuals)
    end
    for i in 1:nworkers(etor)
        etor.params_status[i][1] = PEStatus_Shutdown
    end
    etor
end

function update_done_jobs!(etor::ParallelEvaluator, job_id)
    if job_id > etor.max_done_job
        etor.max_done_job = job_id
    end
    if job_id == etor.max_seq_done_job+1
        # the next sequential job
        etor.max_seq_done_job = job_id
        # see if max_seq_done_job could be further advanced using done jobs
        while etor.max_done_job > etor.max_seq_done_job && first(etor.done_jobs) == etor.max_seq_done_job+1
            etor.max_seq_done_job += 1
            shift!(etor.done_jobs)
        end
    else
        push!(etor.done_jobs, job_id)
    end
end

function process_fitness(etor::ParallelEvaluator, worker_ix::Int)
    update_archive!(etor, job_id, new_fitness)
end

function update_archive!(etor::ParallelEvaluator{F}, job_id::Int, fness::F) where F
    # update the list of done jobs
    #@info "update_archive()"
    candi = pop!(etor.waiting_candidates, job_id)
    etor.unclaimed_candidates[job_id] = candi
    @assert length(etor.unclaimed_candidates) <= 1000_000 # sanity check
    update_done_jobs!(etor, job_id)
    etor.last_fitness = fness
    candi.fitness = archived_fitness(fness, etor.archive)
    etor.num_evals += 1
    #@info "update_archive(): add_candidate()"
    add_candidate!(etor.archive, candi.fitness, candi.params, candi.tag, etor.num_evals)
    return
end

"""
Process all incoming "fitness ready" messages until the evaluator is stopped.
"""
function workers_listener!(etor::ParallelEvaluator{F}) where F
    @info "workers_listener!() started"
    while !is_stopping(etor) || !isempty(etor.waiting_candidates)
        # master critical section
        @inbounds for worker_ix in 1:nworkers(etor)
            #@info "workers_listener!(): checking worker #$worker_ix..."
            #@assert check_worker_running(etor.worker_refs[worker_ix])
            if etor.busy_workers[worker_ix] && (job_id = etor.fitnesses_status[worker_ix][1]) != PEStatus_OK
                @assert (job_id > 0 || is_stopping(etor)) "Worker #$worker_ix bad status: $job_id"
                #@info "worker_listener!(): fitness_evaluated"
                lock(etor.job_assignment)
                @inbounds new_fitness = getfitness(F, etor.shared_fitnesses[worker_ix])
                #@info "worker_listener!($worker_ix): got fitness for job #$job_id"
                etor.fitnesses_status[worker_ix][1] = PEStatus_OK # received
                etor.busy_workers[worker_ix] = false # available again

                unlock(etor.job_assignment)
                Base.release(etor.fitness_slots)

                param_status = etor.params_status[worker_ix][1]
                if param_status == PEStatus_OK # communication in normal state, update the archive
                    update_archive!(etor, job_id, new_fitness)
                elseif param_status < 0 # error/stopping on the master side
                    # remove the candidate
                    delete!(etor.waiting_candidates, job_id)
                end
                #@info "workers_listener!(): yield to other tasks after archive update"
                #yield() # free slots available, switch to the main task
            end
        end
        if length(etor.waiting_candidates) < nworkers(etor)
            if !is_stopping(etor) && isempty(etor.waiting_candidates)
                wait(etor.job_assignment.cond_wait)
            else
                #@info "workers_listener!(): yield to other tasks"
                if !isempty(fitness_done(etor).waitq)
                    # somebody still waiting, notify
                    notify(fitness_done(etor))
                end
                yield() # free slots available, switch to the main task
            end
        end
    end
    @info "workers_listener!() stopped"
end

"""
Asynchronously calculate the fitness of a candidate.
If `force`, existing fitness would be re-evaluated.

Returns -2 if evaluator does not accept jobs because it's shutting down
        -1 if no fitness evaluation was scheduled (`wait=false` and all workers occupied),
        0 if fitness is already evaluated,
        id of fitness evaluation job (check status using `isready()`)
"""
function async_update_fitness(
        etor::ParallelEvaluator{F,FA}, candi::Candidate{FA};
        force::Bool=false, wait::Bool=false) where {F, FA}
    #@info "async_update_fitness(): starting to assign job #$(etor.next_job_id)"
    job_id = etor.next_job_id # tentative job id, but not assigned yet
    if etor.is_stopping
        return -2 # doesn't accept jobs
    elseif !force && !isnafitness(fitness(candi), fitness_scheme(etor.archive))
        return 0 # the candidate has fitness, skip recalculation
    end
    if length(etor.waiting_candidates) >= queue_capacity(etor) && !wait
        #@info "async_update_fitness(): queue is full, skip"
        return -1 # queue full, job not submitted
    end
    #@info "async_update_fitness(): waiting to assign job #$(etor.next_job_id)"
    Base.acquire(etor.fitness_slots)
    #@info "async_update_fitness(): initial slot_state: $(etor.busy_workers), $(etor.fitness_slots.curr_cnt)"
    lock(etor.job_assignment)
    worker_ix = Base.findfirstnot(etor.busy_workers)
    @assert worker_ix !== nothing "Cannot find a worker #$(worker_ix) to put a job to"
    etor.busy_workers[worker_ix] = true # busy
    etor.next_job_id += 1
    copyto!(etor.shared_params[worker_ix], candi.params) # share candidate with the workers
    #@info "async_update_fitness(): assigning job #$job_id to worker #$worker_ix"
    etor.waiting_candidates[job_id] = candi
    #@info "async_update_fitness(): assert fitness status"
    #@assert etor.fitnesses_status[worker_ix][1] == PEStatus_OK
    #@assert etor.params_status[worker_ix][1] == PEStatus_OK
    #@info "async_update_fitness(): flip param status"
    etor.params_status[worker_ix][1] = job_id # announce a message (status = job_id)
    #@info "async_update_fitness(): assigned job #$job_id to worker #$worker_ix"
    #@info "async_update_fitness(): unlock job assignment"
    unlock(etor.job_assignment)
    #@info "async_update_fitness(): yield()"
    #yield() # dispatch the job ASAP, without this it's not getting queued
    return job_id
end

"""
    isready(etor::ParallelEvaluator, fit_job_id::Int)

Check if given asynchronous fitness job calculation is complete.
`fit_job_id` is assigned by `async_update_fitness()`.
"""
function Base.isready(etor::ParallelEvaluator{F,FA}, fit_job_id::Int) where {F, FA}
    fit_job_id > 0 || throw(ArgumentError("Incompatible fitness job Id"))
    pop!(etor.unclaimed_candidates, fit_job_id,
         Candidate{FA}(Individual(), -1, etor.arch_nafitness)) # job was claimed
    return fit_job_id <= etor.max_seq_done_job || in(fit_job_id, etor.done_jobs)
end

"""
    process_completed!(f::Function, etor::ParallelEvaluator)

Processes all completed but not yet claimed candidates.
`f` accepts the completed fitness job Id and corresponding candidate,
returns `true` if the candidate was successfully claimed.
"""
function process_completed!(f::Function, etor::ParallelEvaluator)
    for (job_id, candi) in etor.unclaimed_candidates
        if f(job_id, candi)
            # remove job_id from the waiting list and from the unclaimed list
            #@info "process_completed!($job_id)"
            delete!(etor.unclaimed_candidates, job_id)
        end
    end
    return etor
end

"""
    fitness_done(etor::ParallelEvaluator)

Get the condition that is triggered each time fitness evaluation completes.
"""
@inline fitness_done(etor::ParallelEvaluator) = etor.fitness_slots.cond_wait

"""
    update_fitness!(etor, candidates; [force=false])

Calculate fitness of given `candidates`.
Waits until all fitnesses have been calculated.
`force` specifies whether to re-evaluate fitnesses already stored in `candidates`.
"""
function update_fitness!(etor::ParallelEvaluator{F,FA},
                         candidates::Vector{Candidate{FA}};
                         force::Bool=false) where {F,FA}
    # submit the jobs
    job_ids = sizehint!(BitSet(), length(candidates))
    n_pending = 0
    for candi in candidates
        job_id = async_update_fitness(etor, candi, force=force, wait=true)
        #@info "update_fitness!(): got job id #$job_id"
        if job_id > 0
            n_pending += 1
            push!(job_ids, job_id)
        elseif (job_id < 0) || (force && job_id == 0)
            @warn "fitness calculation rejected"
        end
    end
    # wait until it's done and the evaluator is active
    while n_pending > 0 && !is_stopping(etor) &&
          !(isempty(etor.waiting_candidates) && isempty(etor.unclaimed_candidates))
        #@info "job_ids: $job_ids"
        # pick up the candidates that are for us
        process_completed!(etor) do job_id, candi
            our_job = pop!(job_ids, job_id, 0)>0
            our_job && (n_pending -= 1)
            return our_job
        end
        # wait until another fitness calculation event
        if n_pending > 0 && isempty(etor.unclaimed_candidates)
            #@info "update_fitness!(): wait()"
            wait(fitness_done(etor))
        end
    end
    @assert (n_pending == 0) "Fitnesses not evaluated (#$job_ids)"
    return candidates
end

# WARNING it's not efficient to synchronously calculate single fitness using
# asynchronous `ParallelEvaluator`
function fitness(params::Individual, etor::ParallelEvaluator{F,FA}) where {F, FA}
    candi = Candidate{FA}(params, -1, etor.arch_nafitness)
    job_id = async_update_fitness(etor, candi, wait=true)
    @assert job_id > 0
    while !is_stopping(etor) &&
          !(isempty(etor.waiting_candidates) && isempty(etor.unclaimed_candidates))
        if isready(etor, job_id)
            #@info "fitness(): done"
            return fitness(candi)
        else
            #@info "fitness(): wait()"
            wait(fitness_done(etor))
        end
    end
    error("Fitness not evaluated")
end
