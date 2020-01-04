#= master <-> worker params_status/fitness_status codes =#

const MEStatus_Available = 0 # message received; positive statuses = jobid
const MEStatus_Assigned = -1
const MEStatus_Working = -2
const MEStatus_Success = -3
const MEStatus_Finished = -4
const MEStatus_Error = -5

mutable struct MTEvaluatorWorker{FA}
    status::Int
    jobid::Int
    task::Task
    jobid_chan::Channel{Int}
    candi::Union{Candidate{FA}, Nothing}
    num_evals::Int

    MTEvaluatorWorker{FA}(task::Task, ch::Channel{Int}) where FA =
        new{FA}(MEStatus_Available, 0, task, ch, nothing, 0)
end

isbusy(worker::MTEvaluatorWorker) = worker.status == MEStatus_Working || worker.status == MEStatus_Assigned
isdone(worker::MTEvaluatorWorker) = worker.status == MEStatus_Error || worker.status == MEStatus_Success
isavailable(worker::MTEvaluatorWorker) = worker.status == MEStatus_Available
isrunning(worker::MTEvaluatorWorker) = !istaskdone(worker.task)

function assign!(worker::MTEvaluatorWorker, jobid::Integer, candi::Candidate)
    @assert worker.status == MEStatus_Available
    @assert worker.candi === nothing
    @assert worker.jobid == 0
    worker.status = MEStatus_Assigned
    worker.candi = candi
    put!(worker.jobid_chan, jobid)
    return worker
end

"""
Fitness evaluator that asynchronously distributes calculation
among several worker threads.
"""
mutable struct MultithreadEvaluator{F, FA, FS, P<:OptimizationProblem, A<:Archive} <: AbstractAsyncEvaluator{P}
    problem::P          # optimization problem
    archive::A          # archive where good candidates are automatically stored
    num_evals::Int      # fitness evaluations counter
    last_fitness::F     # last fitness
    arch_nafitness::FA  # NA fitness

    workers::Vector{MTEvaluatorWorker{FA}}

    #new_result::Condition            # new result is available
    #results_lock::ReentrantLock      # results update critical section
    results::Dict{Int, Candidate{FA}} # yet unclaimed candidates with calculated fitness

    done_workers::Channel{Int}
    #busy_workers::Base.Semaphore   # gets acquired when a worker needs to be assigned to a task;
    #                               # used to organize waiting when all workers are busy
    #job_assignment::ReentrantLock  # workers assignment critical section

    done_jobids::SlidingBitset      # ids of completed jobs
    next_jobid::Int                 # ID to assign for the next job

    is_stopping::Bool       # whether the evaluator is in the shutdown sequence
    jobids_pool::Vector{BitSet}     # pool of temporary jobids sets for update_fitness!()

    function MultithreadEvaluator(
        problem::P, archive::A;
        nworkers::Integer = Threads.nthreads() - 1
    ) where {P<:OptimizationProblem, A<:Archive}
        nworkers > 0 || throw(ArgumentError("nworkers must be positive"))
        nworkers < Threads.nthreads() ||
            throw(ArgumentError("nworkers must be less than threads count ($(Threads.nthreads()))"))
        fs = fitness_scheme(problem)
        F = fitness_type(fs)
        FA = fitness_type(archive)

        eval = new{F, FA, typeof(fs), P, A}(
            problem, archive,
            0, nafitness(fs), nafitness(FA),
            Vector{MTEvaluatorWorker{FA}}(),
            #Condition(), ReentrantLock(),
            Dict{Int, Candidate{FA}}(),
            Channel{Int}(10*nworkers), #Base.Semaphore(nworkers), #ReentrantLock(),
            SlidingBitset(), 1,
            false,
            Vector{BitSet}()
        )
        create_workers!(eval, nworkers)

        #finalizer(eval, _shutdown!)
        return eval
    end
end

MultithreadEvaluator(
    problem::OptimizationProblem;
    nworkers::Integer = Threads.nthreads() - 1,
    archiveCapacity::Integer = 10) =
    MultithreadEvaluator(problem, TopListArchive(fitness_scheme(problem), numdims(problem), archiveCapacity),
                         nworkers=nworkers)

num_evals(eval::MultithreadEvaluator) = eval.num_evals

archfitness_type(::Type{<:MultithreadEvaluator{F,FA}}) where {F, FA} = FA
archfitness_type(eval::MultithreadEvaluator) = archfitness_type(typeof(eval))

nworkers(eval::MultithreadEvaluator) = length(eval.workers)
#queue_capacity(eval::MultithreadEvaluator) = nworkers(eval)
#nbusyworkers(eval::MultithreadEvaluator) = eval.busy_workers.curr_cnt

is_stopping(eval::MultithreadEvaluator) = eval.is_stopping

# Create and run the evaluator worker.
# It runs in the separate thread than the master process.
function run_mteval_worker(
    eval::MultithreadEvaluator,
    workerix::Int,
    jobid_chan::Channel{Int}
)
    @debug "Initializing MultithreadEvaluator worker #$workerix at thread=#$(Threads.threadid())"
    ini_jobid = take!(jobid_chan) # get initial jobid=0
    @debug "worker #$workerix: got initial jobid=$ini_jobid"
    @assert ini_jobid == 0
    try
        worker = eval.workers[workerix]
        while !is_stopping(eval)
            #@debug "worker #$workerix: waiting for jobid..."
            jobid = take!(worker.jobid_chan)
            #@debug "worker #$workerix: got jobid=$jobid"
            jobid > 0 || continue
            @assert worker.status == MEStatus_Assigned
            @assert worker.jobid == 0
            @assert worker.candi !== nothing
            #@debug "worker #$workerix: calculating fitness jobid=$jobid"
            worker.jobid = jobid
            worker.status = MEStatus_Working
            eval.last_fitness = candi_fitness = fitness(params(worker.candi), eval.problem)
            worker.candi.fitness = archived_fitness(candi_fitness, eval.archive)
            worker.num_evals += 1
            # clear busy state and notify completion
            worker.status = MEStatus_Success
            #@debug "worker #$workerix: notifying jobid=$jobid done"
            put!(eval.done_workers, workerix)
        end
    catch ex
        #eval.workers[workerix].status = MEStatus_Error
        put!(eval.done_workers, workerix)
        # send -id to notify about an error and to release
        # the master from waiting for worker readiness
        @warn "Exception at MultithreadEvaluator worker #$workerix" exception=ex
        showerror(stderr, ex, catch_backtrace())
        rethrow(ex)
    end
    @debug "worker #$workerix finished"
    eval.workers[workerix].status = MEStatus_Finished
    nothing
end

# HACK (ab)use julia internals to make sure that workers are spawned on different threads
# HACK "inspired" by enq_work() (base/task.jl) and Channel ctor (base/channels.jl)
function MTEvaluatorWorker(eval::MultithreadEvaluator, workerix::Integer, tid::Integer)
    ch = Channel{Int}(1)
    task = Task(() -> run_mteval_worker(eval, workerix, ch))
    task.sticky = false
    bind(ch, task)
    ccall(:jl_set_task_tid, Cvoid, (Any, Cint), task, tid-1)
    push!(Base.Workqueues[tid], task)
    ccall(:jl_wakeup_thread, Cvoid, (Int16,), (tid - 1) % Int16)
    return MTEvaluatorWorker{archfitness_type(eval)}(task, ch)
end

# creates workers and assigns them to different threads
function create_workers!(eval::MultithreadEvaluator, nworkers::Integer)
    @debug "Initializing $nworkers multithread workers..."
    FA = archfitness_type(eval)
    eval.workers = [MTEvaluatorWorker(eval, i, i+1) for i in 1:nworkers]
    @debug "create_workers!(): sending initial jobid=0 to the workers"
    for worker in eval.workers
        put!(worker.jobid_chan, 0)
    end
    @info "MultithreadEvaluator: $nworkers workers ready"
end

# shutdown the evaluator, automatically called when the error occurs
function shutdown!(eval::MultithreadEvaluator)
    @debug "shutdown!(MultithreadEvaluator)"
    eval.is_stopping && error("Cannot shutdown!(MultithreadEvaluator) twice")
    eval.is_stopping = true
    # notify the workers that they should shutdown (each worker should pick exactly one message)
    for worker in eval.workers
        isopen(worker.jobid_chan) && put!(worker.jobid_chan, -1)
    end
    # resume workers listener if it is waiting for the new jobs
    put!(eval.done_workers, 0)
    # make sure all the workers tasks are done
    for (i, worker) in enumerate(eval.workers)
        if isrunning(worker)
            @debug "shutdown!(MultithreadEvaluator): worker #$i is still running, waiting..."
            wait(worker.task)
            @debug "shutdown!(MultithreadEvaluator): worker #$i finished"
        end
    end
    #@assert nbusyworkers(eval) == 0 "Some workers are still busy" upon abnormal termination might not hold
    @info "shutdown!(MultithreadEvaluator): all $(nworkers(eval)) workers stopped"
    @info "shutdown!(MultithreadEvaluator): function evals per worker: $(join([worker.num_evals for worker in eval.workers], ", "))"
end

# Wait until the first incoming fitness completion (or any other) notification from any worker.
# Completed worker is made available again, the candidate is put to the results,
# and its fitness is put to the archive
function wait_workers!(eval::MultithreadEvaluator)
    #@debug "check_workers!(): taking next done worker..."
    workerix = take!(eval.done_workers)
    if workerix <= 0
        @debug "check_workers!(): got worker #$workerix, stopping"
        return false
    end

    worker = eval.workers[workerix]
    if worker.status != MEStatus_Success # communication not in a normal state
        @warn "check_workers!(): worker #$workerix jobid=#$(worker.jobid) status=$(worker.status), shutting down"
        worker.jobid = 0
        worker.candi = nothing
        shutdown!(eval)
        return false
    end

    #@debug "check_workers!(): worker #$workerix jobid=$(worker.jobid) success"
    @assert worker.candi !== nothing
    candi = worker.candi
    #@debug "check_workers!(): store jobid=#$(worker.jobid) result"
    #lock(eval.results_lock)
    push!(eval.done_jobids, worker.jobid)
    eval.results[worker.jobid] = candi
    if length(eval.results) > 10*nworkers(eval)
        @warn "$(length(eval.results)) unclaimed result(s)"
    end
    #unlock(eval.results_lock)
    #@debug "check_workers!(): notify new result"
    #notify(eval.new_result)
    #lock(eval.job_assignment)
    worker.status = MEStatus_Available
    worker.jobid = 0
    worker.candi = nothing
    #unlock(eval.job_assignment)
    #@debug "check_workers!(): add_candidate(archive)"
    eval.num_evals += 1
    add_candidate!(eval.archive, candi.fitness, candi.params, candi.tag, num_evals(eval))
    #@debug "check_workers!(): add_candidate(archive) done"
    return true
end

function async_update_fitness(
    eval::MultithreadEvaluator{F,FA}, candi::Candidate{FA};
    force::Bool=false, wait::Bool=false
) where {F, FA}
    jobid = eval.next_jobid # tentative job id, but not assigned yet, only for logging
    @debug "async_update_fitness(jobid=$jobid?): starting to assign job"
    if is_stopping(eval)
        return -2 # doesn't accept jobs
    elseif !force && !isnafitness(fitness(candi), fitness_scheme(eval.archive))
        @debug "async_update_fitness(jobid=$jobid?): don't recalculate fitness, quit"
        return 0 # the candidate has fitness, skip recalculation
    end
    if !any(isavailable, eval.workers)
        if wait
            #@debug "async_update_fitness(): wait for the completed jobs..."
            wait_workers!(eval) || return -2
            #@debug "async_update_fitness(): done waiting"
        else
            @debug "async_update_fitness(jobid=$jobid?): queue is full, quit"
            return -1 # queue full, job not submitted
        end
    end
    #lock(eval.job_assignment)
    #@debug "worker statuses: $([(worker.status, length(worker.jobid_chan.data)) for worker in eval.workers])"
    workerix = Base.findfirst(isavailable, eval.workers)
    #@debug "async_update_fitness(jobid=$jobid?): assigning job to worker #$workerix"
    if workerix === nothing
        #unlock(eval.job_assignment)
        error("Cannot find a worker to put a job to")
    end
    jobid = eval.next_jobid # now assign job id
    #@debug "async_update_fitness(jobid=$jobid): assigning job to worker #$workerix"
    assign!(eval.workers[workerix], jobid, candi)
    @debug "async_update_fitness(jobid=$jobid): assigned job to worker #$workerix"
    eval.next_jobid += 1
    #unlock(eval.job_assignment)

    #@debug "async_update_fitness(jobid=$jobid): job assigned"
    return jobid
end

is_fitness_ready(eval::MultithreadEvaluator, jobid::Integer) =
    jobid > 0 ? in(eval.done_jobids, jobid) :
    throw(ArgumentError("Invalid job Id=$jobid"))

# Processes all unclaimed candidates with calculated fitnesses (results dict).
# `f` accepts the completed fitness job Id and the corresponding candidate,
# returns `true` if the candidate was successfully claimed by `f` and
# removes the candidate from the unclaimed dict.
function claim_calculated!(f::Function, eval::MultithreadEvaluator)
    nclaimed = 0
    for (jobid, candi) in eval.results
        if f(jobid, candi)
            # remove jobid from the results
            #@debug "claim_calculated!(jobid=#$jobid)"
            delete!(eval.results, jobid)
            nclaimed += 1
        end
    end
    return nclaimed
end

function sync_update_fitness!(f::Any, eval::MultithreadEvaluator, jobids::Any)
    while !is_stopping(eval) && !isempty(jobids) &&
          (any(isrunning, eval.workers) || !isempty(eval.results))
        #@debug "sync_update_fitness!(): jobids=$jobids"
        # pick up the candidates that are for us
        claim_calculated!(eval) do jobid, candi
            our_job = pop!(jobids, jobid, 0)>0
            if our_job # job for one of our candidates
                f !== nothing && f(candi) # externally process the candidate
            end
            return our_job
        end
        # wait if all the workers are busy or all the candidates queued
        if !is_stopping(eval) && !isempty(jobids) && any(isrunning, eval.workers)
            #@debug "sync_update_fitness!(): wait for the completed jobs..."
            wait_workers!(eval)
            #@debug "sync_update_fitness!(): waiting done"
        end
    end
    return nothing
end

function update_fitness!(f::Any, eval::MultithreadEvaluator, candidates::Any;
                         force::Bool=false)
    # submit the jobs
    isempty(candidates) && return candidates
    jobids = isempty(eval.jobids_pool) ? BitSet() : empty!(pop!(eval.jobids_pool))
    (Base.IteratorSize(candidates) === Base.HasLength()) && sizehint!(jobids, length(candidates))
    next = iterate(candidates)
    n_queued = 0
    n_processed = 0
    while ((next !== nothing) || n_queued > 0) && !is_stopping(eval)
        #@debug "update_fitness!(): jobids=$jobids"
        # claim candidates that are for us
        nclaimed = claim_calculated!(eval) do jobid, candi
            our_job = pop!(jobids, jobid, 0)>0
            if our_job # job for one of our candidates
                f !== nothing && f(candi) # externally process the candidate
            end
            return our_job
        end
        n_queued -= nclaimed
        n_processed += nclaimed
        if next !== nothing
            # queue the next candidate
            candi, state = next
            jobid = async_update_fitness(eval, candi, force=force, wait=true)
            #@debug "update_fitness!(): got jobid=#$jobid"
            if jobid > 0
                n_queued += 1
                push!(jobids, jobid)
            elseif (jobid < 0) || (force && jobid == 0)
                @warn "fitness calculation rejected"
            end
            next = iterate(candidates, state)
        end
        # wait if all the workers are busy or all the candidates queued
        if !is_stopping(eval) && n_queued > 0 &&
           (next === nothing || !any(isavailable, eval.workers))
            #@debug "update_fitness!(): wait for the completed jobs..."
            wait_workers!(eval)
            #@debug "update_fitness!(): waiting done"
        end
    end
    push!(eval.jobids_pool, jobids)
    @assert (Base.IteratorSize(typeof(candidates)) === Base.HasLength() ?
             n_processed == length(candidates) :
             n_queued == 0) "Fitnesses not evaluated ($jobids)"
    return candidates
end

# WARNING it's not efficient to synchronously calculate single fitness using
# asynchronous `MultithreadEvaluator`
function fitness(params::Individual, eval::MultithreadEvaluator{F,FA}, tag::Int=0) where {F, FA}
    candi = Candidate{FA}(params, -1, eval.arch_nafitness, nothing, tag)
    jobid = async_update_fitness(eval, candi, wait=true)
    @assert jobid > 0
    @debug "fitness(): is_stopping=$(is_stopping(eval)) has_busy_workers=$(any(isrunning, eval.workers)) has_results=$(!isempty(eval.results))"
    while !is_stopping(eval) && !all(isavailable, eval.workers)
        if is_fitness_ready(eval, jobid)
            break
        else
            #@debug "fitness(): job #$jobid wait for workers"
            wait_workers!(eval)
            #@debug "fitness(): job #$jobid resumed"
        end
        #@debug "fitness(): is_stopping=$(is_stopping(eval)) has_busy_workers=$(any(isrunning, eval.workers)) has_results=$(!isempty(eval.results))"
    end
    if is_fitness_ready(eval, jobid)
        #@debug "fitness(): job #$jobid done"
        pop!(eval.results, jobid) # remove from the results, throws if jobid is not there
        return fitness(candi)
    else
        error("Fitness not evaluated")
    end
end
