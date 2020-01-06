# MultithreadEvaluator worker
# that calculates fitnesses in it's own separate thread
mutable struct MTEvaluatorWorker
    task::Task
    num_evals::Int

    MTEvaluatorWorker(task::Task) = new(task, 0)
end

# check whether the worker task is running
isactive(worker::MTEvaluatorWorker) =
    istaskstarted(worker.task) && !(istaskdone(worker.task) || istaskfailed(worker.task))

# Job object returned by async_update_fitness() method of MultithreadEvaluator
mutable struct MTFitnessEvalJob{FA, NFA, I, S} <: AbstractFitnessEvaluationJob{FA}
    candidates::I   # candidates iterator
    next_candidate::Union{Nothing, Tuple{Candidate{FA}, S}} # next candidate and iter state
    nafitness::NFA

    npending::Threads.Atomic{Int}   # how many candidates are currently being evaluated
    results_lock::Threads.SpinLock
    results::Vector{Candidate{FA}}
    discarded::Bool
end

# all candidates have been iterated over by workers
has_next_candidate(job::MTFitnessEvalJob) = !isnothing(job.next_candidate)
# fitness calculation of all candidates is done
iscomplete(job::MTFitnessEvalJob) = isnothing(job.next_candidate) && (job.npending[] == 0)
Base.isready(job::MTFitnessEvalJob) = iscomplete(job)

# all calculated fitness candidates are claimed by the caller
isclaimed(job::MTFitnessEvalJob) = iscomplete(job) && isempty(job.results)

# wrapper around iterate()
iterate_candidates(candidates, state, nafitness::Nothing) =
    isnothing(state) ? iterate(candidates) : iterate(candidates, state)

# wrapper around iterate() that skips candidate, if their fitness
# doens't match nafitness
function iterate_candidates(candidates, state, nafitness)
    next = isnothing(state) ? iterate(candidates) : iterate(candidates, state)
    while !isnothing(next) && !isequal(next[1].fitness, nafitness)
        next = iterate(candidates, next[2])
    end
    return next
end

# pick the next candidate from the iterator, requires Evaluator's jobs_queue_lock
function take_candidate!(job::MTFitnessEvalJob)
    isnothing(job.next_candidate) && error("No candidates to take")
    candi = job.next_candidate[1]
    # increment before iterating to avoid job being complete for a short while
    Threads.atomic_add!(job.npending, 1)
    # prepare next candidate
    job.next_candidate = iterate_candidates(job.candidates, job.next_candidate[2], job.nafitness)
    return candi
end

# put the candidate to the results, requires job's results_lock
function put_candidate!(job::MTFitnessEvalJob, candi::Candidate)
    job.npending[] > 0 || error("Job doesn't expect candidates to be put back")
    push!(job.results, candi)
    # decrement after adding to avoid job being subject for discarding
    Threads.atomic_sub!(job.npending, 1)
    return job
end

"""
Fitness evaluator that asynchronously distributes calculation
among several worker threads.
"""
mutable struct MultithreadEvaluator{F, FA, FS, P<:OptimizationProblem, A<:Archive} <: AbstractAsyncEvaluator{P}
    problem::P          # optimization problem
    archive::A          # archive where good candidates are automatically stored
    num_evals::Int      # fitness evaluations counter
    num_jobs::Int       # counter of completed claimed jobs
    last_fitness::F     # last fitness
    arch_nafitness::FA  # NA fitness

    is_stopping::Bool       # whether the evaluator is in the shutdown sequence
    workers::Vector{MTEvaluatorWorker}

    jobs_queue::Vector{MTFitnessEvalJob} # jobs for workers to pick up
    jobs_queue_lock::Threads.SpinLock
    iterated_jobs::Vector{MTFitnessEvalJob} # jobs with all candidates iterated (some calculations might be pending), but not yet fully claimed
    iterated_jobs_lock::Threads.SpinLock

    resources_lock::Threads.SpinLock
    atomics_pool::Vector{Threads.Atomic{Int}}
    locks_pool::Vector{Threads.SpinLock}
    results_pool::Vector{Vector{Candidate{FA}}}

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
            0, 0, nafitness(fs), nafitness(FA), false,
            Vector{MTEvaluatorWorker}(),
            Vector{MTFitnessEvalJob}(), Threads.SpinLock(),
            Vector{MTFitnessEvalJob}(), Threads.SpinLock(),
            Threads.SpinLock(),
            Vector{Threads.Atomic{Int}}(),
            Vector{Threads.SpinLock}(),
            Vector{Vector{Candidate{FA}}}(),
        )
        create_workers!(eval, nworkers)

        finalizer(_shutdown!, eval)
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

# FIXME move these method to abstract Evaluator (it would need to support A and FA)
archfitness_type(::Type{<:MultithreadEvaluator{F,FA}}) where {F, FA} = FA
archfitness_type(eval::MultithreadEvaluator) = archfitness_type(typeof(eval))
candidate_type(::Type{T}) where T<:MultithreadEvaluator = Candidate{archfitness_type(T)}
candidate_type(eval::MultithreadEvaluator) = candidate_type(typeof(eval))

nworkers(eval::MultithreadEvaluator) = length(eval.workers)

is_stopping(eval::MultithreadEvaluator) = eval.is_stopping

# what the worker thread should do while waiting for a lock/condition
pause(worker::MTEvaluatorWorker; longer::Bool=false) =
    longer ? yield() : ccall(:jl_cpu_pause, Cvoid, ())

# what the master evaluator thread should do while waiting for a lock/condition
pause(eval::MultithreadEvaluator; longer::Bool=false) =
    longer ? yield() : ccall(:jl_cpu_pause, Cvoid, ());

# Locks SpinLock in a "smart" way: allows yielding to other task from time to
# time, cancels waiting for a lock if the evaluator is being shut down.
# Returns true if the lock was acquired
function smartlock(spinlock::Threads.SpinLock, eval::MultithreadEvaluator,
                   worker::Union{MTEvaluatorWorker, Nothing} = nothing;
                   adaptive::Bool=true)
    locked = false
    threadowner = isnothing(worker) ? eval : worker
    if adaptive
        i = 0
        while !is_stopping(eval) && !(locked = trylock(spinlock))
            pause(threadowner, longer=((i += 1) & 0xFFFF) == 0)
        end
    else
        while !is_stopping(eval) && !(locked = trylock(spinlock))
            pause(threadowner, longer=false)
        end
    end
    return locked
end

# constructs Job object given the candidates iterator.
# Skips candidates with non-natifness, if nafitness !== nothing.
function MTFitnessEvalJob(eval::MultithreadEvaluator, candidates::Any, nafitness::Any)
    next_candidate = iterate_candidates(candidates, nothing, nafitness)
    isnothing(next_candidate) && return nothing

    if smartlock(eval.resources_lock, eval)
        res = MTFitnessEvalJob{archfitness_type(eval), typeof(nafitness),
                               typeof(candidates), typeof(last(next_candidate))}(
            candidates, next_candidate, nafitness,
            isempty(eval.atomics_pool) ? Threads.Atomic{Int}(0) : pop!(eval.atomics_pool),
            isempty(eval.locks_pool) ? Threads.SpinLock() : pop!(eval.locks_pool),
            isempty(eval.results_pool) ? eltype(eval.results_pool)() : pop!(eval.results_pool),
            false)
        unlock(eval.resources_lock)
        return res
    else
        error("Failed to lock MultithreadEvaluator resources")
    end
end

# create MTEvaluatorWorker and assign it to a given tread
# HACK (ab)use julia internals to make sure that workers are spawned on different threads
# HACK "inspired" by enq_work() (base/task.jl) and Channel ctor (base/channels.jl)
function MTEvaluatorWorker(eval::MultithreadEvaluator, workerix::Integer, tid::Integer,
                           readysteadygo::Threads.Atomic{Int})
    task = Task(() -> run_mteval_worker(eval, workerix, readysteadygo))
    task.sticky = true
    ccall(:jl_set_task_tid, Cvoid, (Any, Cint), task, tid-1)
    push!(Base.Workqueues[tid], task)
    ccall(:jl_wakeup_thread, Cvoid, (Int16,), (tid - 1) % Int16)
    return MTEvaluatorWorker(task)
end

# creates workers and assigns them to different threads
function create_workers!(eval::MultithreadEvaluator, nworkers::Integer)
    @debug "Initializing $nworkers multithread workers..."
    readysteadygo = Threads.Atomic{Int}(-nworkers)
    eval.workers = [MTEvaluatorWorker(eval, i, i >= Threads.threadid() ? i+1 : i, readysteadygo)
                    for i in 1:nworkers]
    @debug "Waiting for workers initialization..."
    while !any(w -> istaskfailed(w.task) || istaskdone(w.task), eval.workers) &&
        (readysteadygo[] < 0)
        pause(eval, longer=true)
    end
    if !all(isactive, eval.workers) || readysteadygo[] < 0
        @info "MultithreadEvaluator: workers initialization failed, shutting down"
        shutdown!(eval)
    else
        @info "MultithreadEvaluator: $nworkers workers ready, starting"
        Threads.atomic_add!(readysteadygo, 1)
    end
    return eval.workers
end

# Runs the evaluator worker until the evaluator is shut down.
# The workers look and pick up the first job in the jobs_queue and send it
# to calculate_fitness()
function run_mteval_worker(
    eval::MultithreadEvaluator,
    workerix::Int,
    readysteadygo::Threads.Atomic{Int}
)
    @debug "Initializing MultithreadEvaluator worker #$workerix at thread=#$(Threads.threadid())"
    Threads.atomic_add!(readysteadygo, 1)
    while !is_stopping(eval) && readysteadygo[] <= 0
        yield()
    end
    @debug "MultithreadEvaluator worker #$workerix started"
    try
        worker = eval.workers[workerix]
        i = 0
        nlock_failed = 0
        while !is_stopping(eval)
            queuelock = false
            if isempty(eval.jobs_queue) || islocked(eval.jobs_queue_lock)
                nlock_failed = 0
            else
                if !(queuelock = trylock(eval.jobs_queue_lock))
                    if (nlock_failed += 1) >= 1000000
                        nlock_failed = 0
                        @warn("worker #$workerix: waiting for job for too long (jobs_queue=$(length(eval.jobs_queue)), iterated_jobs=$(length(eval.iterated_jobs)), worker_num_evals=$(worker.num_evals), total_num_evals=$(num_evals(eval)), num_jobs=$(eval.num_jobs))")
                    end
                end
            end
            if queuelock
                nlock_failed = 0
                if !is_stopping(eval) && !isempty(eval.jobs_queue)
                    #@debug "worker #$workerix: got job, calculating"
                    calculate_fitness(first(eval.jobs_queue), eval, workerix)
                    # jobs_queue unlocked by calculate_fitness()
                else
                    @debug "worker #$workerix: locked empty jobs queue ($(length(eval.jobs_queue))) or stopping evaluator ($(is_stopping(eval))), worker_num_evals=$(worker.num_evals), total_num_evals=$(num_evals(eval)), num_jobs=$(eval.num_jobs)"
                    unlock(eval.jobs_queue_lock)
                end
            else
                pause(worker, longer=((i += 1) & 0xFFFF) == 0)
            end
        end
    catch ex
        @warn "Exception at MultithreadEvaluator worker #$workerix" exception=ex
        showerror(stderr, ex, catch_backtrace())
        rethrow(ex)
    end
    @debug "MultithreadEvaluator worker #$workerix finished"
    nothing
end

# Calculates the fitness of the next candidate in the job and puts the
# candidate with the evaluated fitness to job.results.
# Moves the job from eval.jobs_queue to eval.iterated_jobs if it has no
# other candidates.
function calculate_fitness(
    job::MTFitnessEvalJob,
    eval::MultithreadEvaluator,
    workerix::Int
)
    worker = eval.workers[workerix]
    #@debug "worker #$workerix calculate(): taking candidate"
    candi = take_candidate!(job)
    @assert !isnothing(candi)
    no_next = !has_next_candidate(job)
    if no_next
        #@debug "worker #$workerix calculate(): removing iterated job from the jobs_queue"
        first_job = popfirst!(eval.jobs_queue)
        @assert first_job === job # the job should be the first in the running list
    end
    unlock(eval.jobs_queue_lock) # release the jobs queue lock as soon as possible

    if no_next
        # move the job with no further candidates to iterated
        #@debug "worker #$workerix calculate(): move the job to iterated"
        if smartlock(eval.iterated_jobs_lock, eval, worker)
            push!(eval.iterated_jobs, job)
            if length(eval.iterated_jobs) > 100
                @warn "$(length(eval.iterated_jobs)) unclaimed iterated job(s)"
            end
            unlock(eval.iterated_jobs_lock)
        elseif !is_stopping(eval)
            error("Failed to lock iterated_jobs")
        end
    end
    #@debug "worker #$workerix calculate(): calculating fitness..."
    eval.last_fitness = candi_fitness = fitness(params(candi), eval.problem)
    candi.fitness = archived_fitness(candi_fitness, eval.archive)
    worker.num_evals += 1

    #@debug "worker #$workerix calculate(): putting the result back to the job..."
    if smartlock(job.results_lock, eval, worker)
        put_candidate!(job, candi)
        unlock(job.results_lock)
    elseif !is_stopping(eval)
        error("Failed to lock job results")
    end
    return nothing
end

# silent (no IO) shutdown of the evaluator, to please the finalizer()
# set the is_stopping flag so that the workers can quit from their locks
function _shutdown!(eval::MultithreadEvaluator)
    eval.is_stopping = true
end

# shutdown the evaluator, automatically called when the error occurs
function shutdown!(eval::MultithreadEvaluator)
    @debug "shutdown!(MultithreadEvaluator)"
    eval.is_stopping && error("Cannot shutdown!(MultithreadEvaluator) twice")
    eval.is_stopping = true
    # make sure all the workers tasks are done
    for (i, worker) in enumerate(eval.workers)
        if isactive(worker)
            @debug "shutdown!(MultithreadEvaluator): worker #$i is still running, waiting..."
            wait(worker.task)
            @debug "shutdown!(MultithreadEvaluator): worker #$i finished"
        end
    end
    @info "shutdown!(MultithreadEvaluator): all $(nworkers(eval)) workers stopped"
    @info "shutdown!(MultithreadEvaluator): function evals per worker: $(join([worker.num_evals for worker in eval.workers], ", "))"
end

# Discards the job object releasing the reusable resources back to the
# evaluator's resources pool.
# Called from sync_update_fitness()
function discard!(eval::MultithreadEvaluator, job::MTFitnessEvalJob)
    #@debug "discard!(): discarding job..."
    job.discarded && error("Job is already discarded")
    @assert isclaimed(job)
    job.discarded = true

    #@debug "discard!(): removing the job from iterated_jobs"
    if smartlock(eval.iterated_jobs_lock, eval)
        ix = findfirst(j -> j===job, eval.iterated_jobs)
        @assert !isnothing(ix)
        deleteat!(eval.iterated_jobs, ix)
        unlock(eval.iterated_jobs_lock)
    elseif !is_stopping(eval)
        error("Failed to lock iterated_jobs")
    end

    #@debug "discard!(): return the job resources to the pool"
    if smartlock(eval.resources_lock, eval)
        push!(eval.results_pool, empty!(job.results))
        @assert !islocked(job.results_lock)
        @assert job.npending[] == 0
        push!(eval.locks_pool, job.results_lock)
        push!(eval.atomics_pool, job.npending)
        unlock(eval.resources_lock)
    elseif !is_stopping(eval)
        error("Failed to lock resources")
    end
    #@debug "discard!(): done"
    return nothing
end

# Creates MTFitnessEvalJob object and puts it to the jobs_queue
function async_update_fitness!(eval::MultithreadEvaluator, candidates::Any; force::Bool = false)
    #@debug "async_update_fitness!(): creating job..."
    job = MTFitnessEvalJob(eval, candidates, force ? nothing : eval.arch_nafitness)
    isnothing(job) && return job # empty job, do not process further
    #@debug "async_update_fitness!(): adding the job to the queue"
    if smartlock(eval.jobs_queue_lock, eval)
        push!(eval.jobs_queue, job)
        unlock(eval.jobs_queue_lock)
    elseif !is_stopping(eval)
        error("Failed to lock jobs_queue")
    end
    #@debug "async_update_fitness!(): done"
    return job
end

sync_update_fitness(f::Any, job::Nothing, eval::MultithreadEvaluator) = nothing

# Waits until MTFitnessEvalJob is fully processed.
# Constantly monitors job.results for the new candidates and adds them
# to the fitnesses archive.
function sync_update_fitness(f::Any, job::MTFitnessEvalJob, eval::MultithreadEvaluator)
    while !is_stopping(eval) && !isclaimed(job) # continue until all candidates processed
        #@debug "sync_update_fitness(): wait until new results: hasnext=$(has_next_candidate(job)), npending=$(job.npending[]), results=$(length(job.results))"
        i = 0
        while !is_stopping(eval) && isempty(job.results)
            if ((i+=1) & 0xFFFF) == 0
                if all(isactive, eval.workers)
                    pause(eval, longer=true)
                else
                    @warn "MultithreadEvaluator: some workers have stopped, shutting down"
                    shutdown!(eval)
                end
            else
                pause(eval, longer=false)
            end
        end
        is_stopping(eval) && continue
        #@debug "sync_update_fitness(): locking results and resources"
        reslock = false
        if (reslock = smartlock(job.results_lock, eval)) && smartlock(eval.resources_lock, eval)
            #@debug "sync_update_fitness(): swapping results with an empty vector"
            @assert !isempty(job.results)
            results = job.results
            job.results = isempty(eval.results_pool) ?
                similar(results, 0) : pop!(eval.results_pool)
            @assert isempty(job.results)
            unlock(job.results_lock) # the workers can put their results
            reslock = false
            unlock(eval.resources_lock)
            #@debug "sync_update_fitness(): adding $(length(results)) candidate(s) to archive"
            for candi in results
                add_candidate!(eval.archive, candi.fitness, candi.params, candi.tag, eval.num_evals += 1)
                f !== nothing && f(candi)
            end
            #@debug "sync_update_fitness(): returning results vector to the pool"
            if smartlock(eval.resources_lock, eval)
                push!(eval.results_pool, empty!(results))
                unlock(eval.resources_lock)
            elseif !is_stopping(eval)
                error("Failed to lock resources")
            end
        elseif !is_stopping(eval)
            error("Failed to lock resources and job results")
        end
        reslock && unlock(job.results_lock)
    end
    #@debug "sync_update_fitness(): job is done, discarding"
    discard!(eval, job)
    eval.num_jobs += 1

    return job.candidates
end

update_fitness!(f::Any, eval::MultithreadEvaluator, candidates::Any; force::Bool=false) =
    sync_update_fitness(f, async_update_fitness!(eval, candidates, force=force), eval)

update_fitness!(f::Any, eval::MultithreadEvaluator, candidate::Candidate; force::Bool=false) =
    update_fitness!(f, eval, (candidate,), force=force)[1]

# WARNING it's not efficient to synchronously calculate single fitness using
# asynchronous `MultithreadEvaluator`
function fitness(params::Individual, eval::MultithreadEvaluator, tag::Int=0)
    candi = candidate_type(eval)(params, -1, eval.arch_nafitness, nothing, tag)
    update_fitness!(eval, (candi,))
    return fitness(candi)
end
