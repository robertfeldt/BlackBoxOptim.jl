const ParallelPopulationOptimizer_DefaultParameters = @compat Dict{Symbol,Any}(
  :WorkerMethod => :adaptive_de_rand_1_bin_radiuslimited,  # worker population optimization method
  :NWorkers => 2,                   # number of workers
  :Workers => Vector{Int}(),        # IDs of worker processes, takes precedence over NWorkers
  :MigrationSize => 1,              # number of "migrant" individual to sent to the master
  :MigrationPeriod => 100,          # number of worker iterations before sending migrants
  :ArchiveCapacity => 10,           # ParallelPseudoEvaluator archive capacity
  :ToWorkerChannelCapacity => 100,  # how many unread messages the master->worker channel can store
  :FromWorkersChannelCapacity => 1000 # how many unread messages the workers->master channel can store
)

# metrics from worker optimizer
type WorkerMetrics
  num_evals::Int           # number of function evals
  num_steps::Int           # number of steps
  num_better::Int          # number of steps that improved best fitness
  num_migrated::Int        # number of migrants worker received
  num_better_migrated::Int # number of migrants accepted (because fitness improved)

  WorkerMetrics() = new(0, 0, 0, 0, 0)
end

function Base.copy!(a::WorkerMetrics, b::WorkerMetrics)
  a.num_evals = b.num_evals
  a.num_steps = b.num_steps
  a.num_better = b.num_better
  a.num_migrated = b.num_migrated
  a.num_better_migrated = b.num_better_migrated
  return a
end

function Base.show(io::IO, m::WorkerMetrics)
    println(io, "evals=", m.num_evals, " steps=", m.num_steps,
                " improvements=", m.num_better,
                " migrated=", m.num_migrated,
                " (improved=", m.num_better_migrated, ")")
end

# fake evaluator for ParallelPopulationOptimizer
# it doesn't evaluate itself, but stores some
# metrics from the workers evaluators
type ParallelPseudoEvaluator{F, P<:OptimizationProblem} <: Evaluator{P}
  problem::P
  archive::TopListArchive{F}
  workers_metrics::Vector{WorkerMetrics}  # function evals per worker etc
  last_fitness::F
end

num_evals(ppe::ParallelPseudoEvaluator) = mapreduce(x -> x.num_evals, +, 0, ppe.workers_metrics)
num_better(ppe::ParallelPseudoEvaluator) = mapreduce(x -> x.num_better, +, 0, ppe.workers_metrics)

function ParallelPseudoEvaluator{P<:OptimizationProblem}(
    problem::P, nworkers::Int;
    ArchiveCapacity::Int = 10)
    ParallelPseudoEvaluator{fitness_type(problem), P}(
        problem,
        TopListArchive(fitness_scheme(problem), numdims(problem), ArchiveCapacity),
        WorkerMetrics[WorkerMetrics() for i in 1:nworkers],
        nafitness(fitness_scheme(problem)))
end

const FINAL_CANDIDATE = -12345 # special terminating candidate with worker index sent by master
const ERROR_CANDIDATE = -67890 # special error candidate with worker index send by worker

# message with the candidate passed between the workers and the master
immutable CandidateMessage{F}
  worker::Int               # origin of candidate
  metrics::WorkerMetrics    # current origin worker metrics
  candi::Candidate{F}
end

typealias WorkerChannel{F} Channel{CandidateMessage{F}}
typealias WorkerChannelRef{F} RemoteRef{WorkerChannel{F}}

# Parallel population optimizer
# starts nworkers parallel population optimizers.
# At regular interval, the workers send the master process their random population members
# and the master redirects them to the other workers
type ParallelPopulationOptimizer{F, P<:OptimizationProblem} <: SteppingOptimizer
  worker_procs::Vector{Int}                         # IDs of worker processes
  optimizer_generator::Function                     # generates optimizer on worker process
  migrationPeriod::Int                              # worker optimizer iterations between sending the migrations
  migrationSize::Int                                # how many individuals to send at once
  final_fitnesses::Vector{RemoteRef{Channel{Any}}}  # references to the @spawnat ID run_worker()
  from_workers::WorkerChannelRef{F}                 # inbound channel of candidates from all workers
  to_workers::Vector{WorkerChannelRef{F}}           # outgoing channels to each worker
  is_started::RemoteRef{Channel{Bool}}              # flag that all workers have started
  evaluator::ParallelPseudoEvaluator{F, P}          # aggregates workers optimization states
  population::Any # FIXME use base abstract type for population, when it would be available
end

nworkers(ppopt::ParallelPopulationOptimizer) = length(ppopt.worker_procs)
population(ppopt::ParallelPopulationOptimizer) = ppopt.population

# read worker's message, stores the worker metrics and updates best fitness using
function store!{F}(ppe::ParallelPseudoEvaluator{F}, msg::CandidateMessage{F})
  copy!(ppe.workers_metrics[msg.worker], msg.metrics)
  if !isnafitness(msg.candi.fitness, fitness_scheme(ppe)) # store only the candidates with the known fitness
    add_candidate!(ppe.archive, msg.candi.fitness, msg.candi.params, num_evals(ppe))
  end
end

# check that worker is stil running.
# If running, its RemoteRefs should not be ready,
# but if there was exception in the worker,
# it would be thrown into the main thread
function check_worker_running{T}(worker::RemoteRef{T})
  if isready(worker)
    fetch(worker) # fetch the worker, this should trigger an exception
    # no exception, but the worker should not be ready
    error("Worker #? has finished before the master shutdown")
  end
  return true
end

# outer parallel population optimizer constructor that
# also spawns worker tasks
function ParallelPopulationOptimizer{P<:OptimizationProblem}(
    problem::P, optimizer_generator::Function; NWorkers::Int = 1, Workers::Vector{Int} = Vector{Int}(),
    MigrationSize::Int = 1, MigrationPeriod::Int = 100,
    ArchiveCapacity::Int = 10,
    ToWorkerChannelCapacity::Int = 1000,
    FromWorkersChannelCapacity::Int = (isempty(Workers) ? NWorkers : length(Workers)) * ToWorkerChannelCapacity)
  info("Constructing parallel optimizer...")
  F = fitness_type(problem)
  if isempty(Workers)
    # take the first NWorkers workers
    Workers = workers()[1:NWorkers]
  end
  if length(Workers) <= 1
    throw(ArgumentError("ParallelPopulationOptimizer requires at least 2 worker processes"))
  end
  ParallelPopulationOptimizer{F, P}(Workers, optimizer_generator, MigrationPeriod, MigrationSize,
       Vector{RemoteRef{Channel{Any}}}(length(Workers)),
       RemoteRef(() -> WorkerChannel{F}(FromWorkersChannelCapacity)),
       WorkerChannelRef{F}[RemoteRef(() -> WorkerChannel{F}(ToWorkerChannelCapacity), id) for id in Workers],
       RemoteRef(() -> Channel{Bool}(1)),
      ParallelPseudoEvaluator(problem, length(Workers);
                              ArchiveCapacity = ArchiveCapacity), nothing)
end

function setup!(ppopt::ParallelPopulationOptimizer)
  info("Initializing parallel workers...")
  workers_ready = RemoteRef(() -> Channel{Int}(nworkers(ppopt))) # FIXME do we need to wait for the worker?
  @assert !isready(ppopt.is_started)
  # spawn suboptimizers
  for i in eachindex(ppopt.worker_procs)
    procid = ppopt.worker_procs[i]
    info("  Spawning worker #$i at process #$procid...");
    ppopt.final_fitnesses[i] = @spawnat procid run_worker(i,
                           workers_ready, ppopt.is_started,
                           problem(ppopt.evaluator), ppopt.optimizer_generator,
                           ppopt.from_workers, ppopt.to_workers[i],
                           ppopt.migrationSize, ppopt.migrationPeriod)
  end
  # wait until all the workers are started
  info("Waiting for the workers to be ready...")
  # FIXME is it required?
  nready = 0
  while nready < nworkers(ppopt)
    map(check_worker_running, ppopt.final_fitnesses)
    worker_id = take!(workers_ready)
    if worker_id < 0
      # worker failed to initialize, reading its task would throw an exception
      check_worker_running(ppopt.final_fitnesses[-worker_id])
      error("Exception in the worker, but all workers still running")
    end
    info("  Worker #$worker_id is ready")
    nready += 1
  end
  info("All workers ready")
  return ppopt
end

function parallel_population_optimizer(problem::OptimizationProblem, parameters::Parameters)
  param_dict = convert(ParamsDict, parameters) # FIXME convert to dict to avoid serialization problems of DictChain
  params = chain(ParallelPopulationOptimizer_DefaultParameters, parameters)
  worker_method = params[:WorkerMethod]
  info("Using $worker_method as worker method for parallel optimization")
  optimizer_func = ValidMethods[worker_method]

  ParallelPopulationOptimizer(problem, function (id::Int) problem, optimizer_func(problem, param_dict) end,
                              NWorkers = params[:NWorkers], Workers = params[:Workers],
                              MigrationSize = params[:MigrationSize],
                              MigrationPeriod = params[:MigrationPeriod],
                              ArchiveCapacity = params[:ArchiveCapacity],
                              ToWorkerChannelCapacity = params[:ToWorkerChannelCapacity],
                              FromWorkersChannelCapacity = params[:FromWorkersChannelCapacity])
end

# redirects candidate to another worker
function redirect{F}(ppopt::ParallelPopulationOptimizer{F}, msg::CandidateMessage{F})
  # redirect to the other parallel task
  #println("redirecting from $(msg.worker)")
  recv_ix = sample(1:(length(ppopt.to_workers)-1))
  if recv_ix >= msg.worker # index is the origin worker
    recv_ix += 1
  end
  @assert (recv_ix != msg.worker) && (1 <= recv_ix <= nworkers(ppopt))
  msg.candi.op = NO_GEN_OP # reset operation and tag to avoid calling adjust!() out-of-context
  msg.candi.tag = 0
  put!(ppopt.to_workers[recv_ix], msg)
  #println("redirecting done")
end

function step!{F}(ppopt::ParallelPopulationOptimizer{F})
  #println("main#: n_evals=$(num_evals(ppopt.evaluator))")
  if !isready(ppopt.is_started) put!(ppopt.is_started, true) end # if it's the first iteration
  last_better = num_better(ppopt.evaluator)
  candidate = take!(ppopt.from_workers)::CandidateMessage{F}
  if candidate.worker == ERROR_CANDIDATE
    map(check_worker_running, ppopt.final_fitnesses)
    error("Exception in the worker, but all workers still running")
  end
  #println("candidate=$candidate")
  store!(ppopt.evaluator, candidate)
  redirect(ppopt, candidate)
  return num_better(ppopt.evaluator) - last_better
end

# finalize the master: wait for the workers shutdown,
# get their best candidates
function finalize!{F}(ppopt::ParallelPopulationOptimizer{F})
  info("Finalizing parallel optimizer...")
  # send special terminating candidate
  for to_worker in ppopt.to_workers
    put!(to_worker, CandidateMessage{F}(FINAL_CANDIDATE, WorkerMetrics(),
                                        Candidate{F}(Individual())))
  end
  # wait until all threads finish
  # the last candidates being sent are the best in the population
  info("Waiting for the workers to finish...")
  for i in eachindex(ppopt.final_fitnesses)
    msg, pop = fetch(ppopt.final_fitnesses[i])
    msg::CandidateMessage{F}
    @assert msg.worker == i
    store!(ppopt.evaluator, msg) # store the best candidate
    if ppopt.population === nothing
      ppopt.population = deepcopy(pop)
    elseif isa(ppopt.population, PopulationMatrix)
      # FIXME switch by populaton type is ugly, should go aways when
      # there would be single Population base class
      ppopt.population = hcat(ppopt.population, pop)
    else
      append!(ppopt.population, pop)
    end
    info("Worker #$(msg.worker) finished")
  end
  # clean is_started channel
  while isready(ppopt.is_started)
    take!(ppopt.is_started)
  end
  info("Parallel optimizer finished. Metrics per worker: $(ppopt.evaluator.workers_metrics)")
end

# wraps the worker's population optimizer
# and communicates with the master
type PopulationOptimizerWrapper{F,O<:PopulationOptimizer,E<:Evaluator}
  id::Int                                             # worker's Id
  optimizer::O
  evaluator::E
  to_master::WorkerChannelRef{F}   # outgoing candidates
  from_master::WorkerChannelRef{F} # incoming candidates
  migrationSize::Int     # size of the migrating group
  migrationPeriod::Int   # number of iterations between the migrations

  metrics::WorkerMetrics # current metrics
  is_stopping::Bool      # if the optimizer is being shut down
  can_receive::Condition # condition recv_task waits for
  can_run::Condition     # condition run!() task waits for
  recv_task::Task        # task that continuously executes recv_immigrants()

  function Base.call{F,O,E}(::Type{PopulationOptimizerWrapper},
      id::Int, optimizer::O, evaluator::E,
      to_master::WorkerChannelRef{F}, from_master::WorkerChannelRef{F},
      migrationSize::Int = 1, migrationPeriod::Int = 100)
    res = new{F,O,E}(id, optimizer, evaluator,
                     to_master, from_master,
                     migrationSize, migrationPeriod,
                     WorkerMetrics(), false, Condition(), Condition())
    # "background" migrants receiver task
    res.recv_task = @schedule while !res.is_stopping
      wait(res.can_receive)
      recv_immigrants!(res)
      notify(res.can_run)
    end
    return res
  end
end

function send_emigrants{F}(wrapper::PopulationOptimizerWrapper{F})
  pop = population(wrapper.optimizer)
  # prepare the group of emigrants
  migrant_ixs = sample(1:popsize(pop), wrapper.migrationSize, replace=false)
  for migrant_ix in migrant_ixs
    migrant = acquire_candi(pop, migrant_ix)
    # send them outward
    wrapper.metrics.num_evals = num_evals(wrapper.evaluator)
    put!(wrapper.to_master, CandidateMessage{F}(wrapper.id, wrapper.metrics, migrant))
    # FIXME check that we the reuse of candidate does not affect
    # the migrants while they wait to be sent
    # don't release the candidate back to ensure sent data is not affected
    # release_candi(pop, migrant)
  end
end

# receive migrants (called from "background" task)
function recv_immigrants!{F}(wrapper::PopulationOptimizerWrapper{F})
  pop = population(wrapper.optimizer)
  # receive all the immigrants in the channel, if any; do not wait if it's empty
  while isready(wrapper.from_master)
    msg = take!(wrapper.from_master)::CandidateMessage{F}
    if msg.worker == FINAL_CANDIDATE # special index sent by master to indicate termination
      wrapper.is_stopping = true
      break
    end
    # assign migrants to random population indices
    migrant_ix = sample(1:popsize(pop))
    candidates = sizehint!(Vector{candidate_type(pop)}(), 2)
    push!(candidates, acquire_candi(pop, migrant_ix))
    push!(candidates, acquire_candi(pop, msg.candi))
    candidates[end].index = migrant_ix # override the incoming index
    rank_by_fitness!(wrapper.evaluator, candidates)
    wrapper.metrics.num_migrated += 1
    wrapper.metrics.num_better_migrated += tell!(wrapper.optimizer, candidates)
  end
end

# run the wrapper (called in the "main" task)
function run!(wrapper::PopulationOptimizerWrapper)
  while !wrapper.is_stopping
    if istaskdone(wrapper.recv_task)
      error("recv_task has completed prematurely")
    end
    wrapper.metrics.num_steps += 1
    #println("$(wrapper.metrics.num_steps)-th iteration")
    if wrapper.metrics.num_steps % wrapper.migrationPeriod == 0
      send_emigrants(wrapper) # before we started processing
      notify(wrapper.can_receive) # switch to migrants receiving task
      wait(wrapper.can_run)
    end
    # normal ask/tell sequence
    candidates = ask(wrapper.optimizer)
    rank_by_fitness!(wrapper.evaluator, candidates)
    wrapper.metrics.num_better += tell!(wrapper.optimizer, candidates)
  end
  finalize!(wrapper.optimizer, wrapper.evaluator)
end

# returns the candidate message with final metrics and the best candidate
function final_fitness{F}(wrapper::PopulationOptimizerWrapper{F})
  @assert wrapper.is_stopping
  # get the best candidate in the population
  pop = population(wrapper.optimizer)
  best_candi = acquire_candi(pop)
  copy!(best_candi.params, best_candidate(wrapper.evaluator.archive))
  best_candi.fitness = best_fitness(wrapper.evaluator.archive)
  best_candi.index = -1 # we don't know it
  best_candi.tag = 0
  wrapper.metrics.num_evals = num_evals(wrapper.evaluator)
  # return the best candidate and the whole population
  CandidateMessage{F}(wrapper.id, wrapper.metrics, best_candi), pop
end

# Function that the master process spawns at each worker process.
# Creates and run the worker wrapper
function run_worker{F}(id::Int,
                    worker_ready::RemoteRef{Channel{Int}},
                    is_started::RemoteRef{Channel{Bool}},
                    problem::OptimizationProblem,
                    optimizer_generator::Function,
                    to_master::WorkerChannelRef{F},
                    from_master::WorkerChannelRef{F},
                    migrationSize, migrationPeriod)
  info("Initializing parallel optimization worker #$id at task=$(myid())")
  wrapper = nothing
  try
    problem, opt = optimizer_generator(id)
    wrapper = PopulationOptimizerWrapper(id,
                opt,
                ProblemEvaluator(problem),
                to_master, from_master,
                migrationSize, migrationPeriod)
  catch e
    # send -id to notify about an error and to release
    # the master from waiting for worker readiness
    put!(worker_ready, -id)
    rethrow(e)
  end
  # create immigrants receiving tasks=#
  put!(worker_ready, id)
  info("Waiting for the master start signal...")
  fetch(is_started) # wait until the master is started
  info("Running worker #$id")
  try
    run!(wrapper)
  catch e
    # send error candidate to notify about an error and to release
    # the master from waiting for worker messages
    put!(wrapper.to_master,
         CandidateMessage{F}(ERROR_CANDIDATE, WorkerMetrics(),
                             Candidate{F}(Individual())))
    rethrow(e)
  end
  info("Worker #$id stopped")
  final_fitness(wrapper) # return the best fitness
end
