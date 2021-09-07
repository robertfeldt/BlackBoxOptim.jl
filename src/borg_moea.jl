"""
Borg MOEA algorithm.

Based on Hadka & Reed, "Borg: An Auto-Adaptive Many-Objective Evolutionary Computing
Framework", Evol. Comp. 2013
"""
mutable struct BorgMOEA{FS<:FitnessScheme,V<:Evaluator,P<:Population,M<:GeneticOperator,E<:EmbeddingOperator} <: SteppingOptimizer
    evaluator::V
    population::P     # candidates, NOTE index one is always reserved for one archive parent for recombination
    rand_check_order::Vector{Int} # random population check order

    n_restarts::Int
    n_steps::Int
    last_restart_check::Int
    last_restart::Int
    last_wrecombinate_update::Int

    τ::Float64        # tournament size, fraction of the population
    γ::Float64        # recommended population-to-archive ratio
    γ_δ::Float64      # the maximum allowed deviation of the population-to-archive ratio from γ
    min_popsize::Int  # the minimal population size

    recombinate_distr::Categorical # weights of the recombination operators

    # method parameters
    θ::Float64        # restart-discounting coefficient for recombination operator weights
    ζ::Float64        # dampening coefficient for recombination operator weights
    wrecombinate_update_period::Int
    restart_check_period::Int
    max_steps_without_ϵ_progress::Int

    recombinate::Vector{CrossoverOperator} # recombination operators
    recomb_fit_job::Union{AbstractFitnessEvaluationJob, Nothing} # job for fitness calculation of recombined candidates (for AsyncEval version)

    # Set of operators that together define a specific DE strategy.
    select::TournamentSelector{HatCompare{FS}}         # random individuals selector
    modify::M         # operator to mutate frontier element during restarts
    embed::E          # embedding operator

    function BorgMOEA(
        problem::O,
        pop::P, recombinate::Vector{CrossoverOperator},
        modify::M = M(), embed::E = E(), params = EMPTY_PARAMS) where
            {O<:OptimizationProblem, P<:Population,
             M<:GeneticOperator, E<:EmbeddingOperator}
        # NOTE if ϵ-dominance is used, params[:ϵ] has the priority
        fit_scheme = fitness_scheme(problem)
        isa(fit_scheme, TupleFitnessScheme) || throw(ArgumentError("BorgMOEA can only solve problems with `TupleFitnessScheme`"))
        !isempty(recombinate) || throw(ArgumentError("No recombinate operators specified"))
        fit_scheme = EpsBoxDominanceFitnessScheme(fit_scheme, params[:ϵ])
        archive = EpsBoxArchive(fit_scheme, params)
        evaluator = make_evaluator(problem, archive, params)
        new{typeof(fit_scheme),typeof(evaluator),P,M,E}(evaluator, pop, Vector{Int}(), 0, 0, 0, 0, 0,
                params[:τ], params[:γ], params[:γ_δ], params[:PopulationSize],
                Categorical(ones(length(recombinate))/length(recombinate)),
                params[:θ], params[:ζ], params[:OperatorsUpdatePeriod], params[:RestartCheckPeriod],
                params[:MaxStepsWithoutEpsProgress],
                recombinate, nothing,
                TournamentSelector(fit_scheme, ceil(Int, params[:τ]*popsize(pop))), modify, embed)
    end
end

const BorgMOEA_DefaultParameters = chain(EpsBoxArchive_DefaultParameters, ParamsDict(
    :ϵ => 0.1,        # size of the ϵ-box
    :τ => 0.02,       # selection ratio, fraction of population to use for tournament
    :γ => 4.0,        # recommended population-to-archive ratio
    :γ_δ => 0.25,     # the maximum allowed deviation of the population-to-archive ratio from γ
    :θ => 0.9,        # restart-discounting coefficient for recombination operator weights
    :ζ => 1.0,        # dampening coefficient for recombination operator weights
    :RestartCheckPeriod => 1000,
    :OperatorsUpdatePeriod => 100,
    :MaxStepsWithoutEpsProgress => 100
))

function borg_moea(problem::OptimizationProblem, options::Parameters = EMPTY_PARAMS)
    opts = chain(BorgMOEA_DefaultParameters, options)
    fs = fitness_scheme(problem)
    N = numobjectives(fs)
    F = fitness_eltype(fs)
    pop = population(problem, opts, nafitness(IndexedTupleFitness{N,F}), ntransient=1)
    BorgMOEA(problem, pop, CrossoverOperator[DiffEvoRandBin1(chain(DE_DefaultOptions, options)),
                                           SimulatedBinaryCrossover(chain(SBX_DefaultOptions, options)),
                                           SimplexCrossover{3}(chain(SPX_DefaultOptions, options)),
                                           ParentCentricCrossover{2}(chain(PCX_DefaultOptions, options)),
                                           ParentCentricCrossover{3}(chain(PCX_DefaultOptions, options)),
                                           UnimodalNormalDistributionCrossover{2}(chain(UNDX_DefaultOptions, options)),
                                           UnimodalNormalDistributionCrossover{3}(chain(UNDX_DefaultOptions, options))],
            FixedGeneticOperatorsMixture(GeneticOperator[
                                            MutationClock(PolynomialMutation(search_space(problem), chain(PM_DefaultOptions, options)), 1/numdims(problem)),
                                            MutationClock(UniformMutation(search_space(problem)), 1/numdims(problem))], [0.75, 0.25]),
            RandomBound(search_space(problem)), opts)
end

archive(alg::BorgMOEA) = alg.evaluator.archive

set_candidate!(o::BorgMOEA, x0) = set_candidate!(o.population, x0)

set_multi_candidate!(o::BorgMOEA, x0_list) = set_multi_candidate!(o.population, x0_list)

# Take one step of Borg MOEA.
function step!(alg::BorgMOEA)
    if alg.n_steps == 0
        # make sure fitness is calculated for every population member when starting,
        # it guarantees that no Pareto set elements would be lost when continuing
        # optimization and using the previous pareto_frontier() as the starting population
        update_population_fitness!(alg)
    end

    alg.n_steps += 1
    if alg.n_steps >= alg.last_restart_check + alg.restart_check_period
        alg.last_restart_check = alg.n_steps
        # check for restarting conditions
        if (!isempty(archive(alg)) &&
            (abs(popsize(alg.population) - alg.γ * length(archive(alg))) >= alg.γ_δ * length(archive(alg)))) ||
            noprogress_streak(archive(alg), since_restart=true) >=  alg.max_steps_without_ϵ_progress
            restart!(alg)
        end
    end
    if alg.n_steps >= alg.last_wrecombinate_update + alg.wrecombinate_update_period
        update_recombination_weights!(alg)
    end

    prepare_recombination(alg)
    # Select the operators to apply based on their probabilities
    recomb_op_ix = rand(alg.recombinate_distr)
    recombine_individuals!(alg, recomb_op_ix, alg.recombinate[recomb_op_ix])
    return alg
end

# "kernel function" of step() that would specialize to given xover operator type
function recombine_individuals!(alg::BorgMOEA, recomb_op_ix::Int, recomb_op::CrossoverOperator)
    # select parents for recombination
    n_parents = numparents(recomb_op)
    # parent indices
    parent_indices = select(alg.select, alg.population,
                            isempty(archive(alg)) ? n_parents : n_parents-1)
    if !isempty(archive(alg))
        # get one parent from the archive and copy it to the fitrst transient member
        arch_ix = transient_range(alg.population)[1]
        alg.population[arch_ix] = rand_front_elem(archive(alg))
        push!(parent_indices, arch_ix)
    end
    # Crossover parents and target
    children = acquire_candis(alg.population, numchildren(recomb_op))
    apply!(recomb_op, Individual[child.params for child in children],
           zeros(Int, length(children)), alg.population, parent_indices)
    for child in children
        apply!(alg.embed, child.params, alg.population, parent_indices[1])
        reset_fitness!(child, alg.population)
        child.extra = recomb_op
        child.tag = recomb_op_ix
    end
    postprocess_recombined!(alg, children)
end

prepare_recombination(alg::BorgMOEA) = nothing # do nothing

# AsyncEvaluator version -- process previously submitted candidates with the completed fitness
function prepare_recombination(alg::BorgMOEA{<:FitnessScheme, <:AbstractAsyncEvaluator})
    if !isnothing(alg.recomb_fit_job)
        sync_update_fitness(alg.recomb_fit_job, alg.evaluator) do candi
            process_candidate!(alg, candi)
            return true
        end
        @assert isready(alg.recomb_fit_job)
        alg.recomb_fit_job = nothing
    end
    return nothing
end

function postprocess_recombined!(alg::BorgMOEA, candidates::Any)
    update_fitness!(alg.evaluator, candidates, force=true) # implicitly updates the archive
    process_candidate!.(Ref(alg), candidates)
    return candidates
end

# AsyncEvaluator version, just submit to fitness calculation, nothing else
# if the queue is full, waits until some jobs are processed -- that established the
# balance between recombining and fitness evaluation
function postprocess_recombined!(alg::BorgMOEA{<:FitnessScheme, <:AbstractAsyncEvaluator}, candidates::Any)
    @assert isnothing(alg.recomb_fit_job)
    alg.recomb_fit_job = async_update_fitness!(alg.evaluator, candidates, force=true)
    @assert !isnothing(alg.recomb_fit_job)
    return candidates
end

function process_candidate!(alg::BorgMOEA, candi::Candidate)
    ifitness = fitness(candi)
    # test the population
    hat_comp = HatCompare(fitness_scheme(archive(alg)))
    popsz = popsize(alg.population)
    if length(alg.rand_check_order) != popsz
        # initialize the random check order
        alg.rand_check_order = randperm(popsz)
    end
    comp = 0
    candi.index = 0
    # iterate through the population in a random way
    n_checks = rand(min(popsz, 2*alg.select.size):popsz)
    @inbounds for i in 1:n_checks
        # use "modern" Fisher-Yates shuffle to gen random population index
        j = rand(i:popsz)
        ix = alg.rand_check_order[j] # the next random individual
        if j > i
            alg.rand_check_order[j] = alg.rand_check_order[i]
            alg.rand_check_order[i] = ix
        end
        cur_comp = hat_comp(ifitness, fitness(alg.population, ix))[1]
        if cur_comp > 0 # new candidate does not dominate
            comp = cur_comp
            break
        elseif cur_comp < 0 # replace the first dominated
            comp = cur_comp
            candi.index = ix
            # FIXME the population check is stopped when the first candidate dominated
            # by the `child` is found, but since the population might already contain candidates
            # dominated by others, it could be that the `child` is also dominated
            # In Borg paper they do not discuss this situation in detail -- whether the search
            # should continue
            break
        end
    end
    if comp == 0 # non-dominated among the checked ones, replace random checked individual
        candi.index = alg.rand_check_order[rand(1:n_checks)]
    end
    if candi.index > 0
        accept_candi!(alg.population, candi)
    else
        release_candi(alg.population, candi)
    end
    return nothing
end

# trace current optimization state,
# Called by OptRunController trace_progress()
function trace_state(io::IO, alg::BorgMOEA, mode::Symbol)
    println(io, "pop.size=", popsize(alg.population),
                " arch.size=", length(archive(alg)),
                " n.restarts=", alg.n_restarts)
    if mode == :verbose
        # output recombination operator rates
        println(io, "P(recombine):")
        for i in eachindex(alg.recombinate)
            println(io, "  #$i(", alg.recombinate[i], ")=",
                    @sprintf("%.4f",alg.recombinate_distr.p[i]))
        end
    end
end

update_population_fitness!(alg::BorgMOEA) =
    update_fitness!(alg.evaluator,
                    PopulationCandidatesIterator(alg.population,
                                                 alg.population.nafitness)) do candi
        alg.population.fitness[candi.index] = candi.fitness
        release_candi(alg.population, candi)
    end

"""
Update recombination operator probabilities based on the archive tag counts.
"""
function update_recombination_weights!(alg::BorgMOEA)
    op_counts = tagcounts(archive(alg), alg.θ)
    adj_op_counts = sum(values(op_counts)) + length(alg.recombinate)*alg.ζ
    alg.recombinate_distr = Categorical([(get(op_counts, i, 0)+alg.ζ)/adj_op_counts
                                        for i in eachindex(alg.recombinate)])
    alg.last_wrecombinate_update = alg.n_steps
    return alg
end

function acquire_mutant(alg::BorgMOEA, ix::Int, last_nonmutant::Int)
    # select one of the existing individuals (must be from the frontier)
    mutant = acquire_candi(alg.population, rand(1:last_nonmutant))
    mutant.tag = 0 # it's not generated by the crossover
    mutant.index = ix
    reset_fitness!(mutant, alg.population)
    apply!(alg.modify, mutant.params, ix)
    # project using the archived individual as the reference
    apply!(alg.embed, mutant.params, alg.population, rand(1:last_nonmutant))
    return mutant
end

"""
Iterates the mutants.
"""
struct BorgMutantsIterator{P<:PopulationWithFitness, A<:BorgMOEA} <: AbstractPopulationCandidatesIterator{P}
    alg::A
    last_nonmutant::Int

    BorgMutantsIterator(alg::A, last_nonmutant::Int) where A<:BorgMOEA =
        new{typeof(alg.population), A}(alg, last_nonmutant)
end

Base.IteratorSize(::Type{<:BorgMutantsIterator}) = Base.HasLength()
Base.length(it::BorgMutantsIterator) = popsize(it.alg.population) - it.last_nonmutant

Base.iterate(it::BorgMutantsIterator, ix::Integer = it.last_nonmutant) =
    ix < popsize(it.alg.population) ? (acquire_mutant(it.alg, ix+1, it.last_nonmutant), ix+1) : nothing

populate_by_mutants(alg::BorgMOEA, last_nonmutant::Integer) =
    update_fitness!(alg.evaluator, BorgMutantsIterator(alg, last_nonmutant)) do mutant
        accept_candi!(alg.population, mutant)
    end

"""
Restart Borg MOEA.

Resize and refills the population from the archive.
"""
function restart!(alg::BorgMOEA)
    notify!(archive(alg), :restart)
    narchived = length(archive(alg))
    new_popsize = max(alg.min_popsize, ceil(Int, alg.γ * narchived))
    # fill populations with the solutions from the archive
    resize!(alg.population, new_popsize)
    ndelegates = min(narchived, new_popsize)
    # get the indexes of the frontier elements to put into the population
    delegate_ixs = sort!(sample(1:narchived, ndelegates, replace=false))
    delegate_pos = 1 # position of the current delegate index
    for (candi_ix, candi) in enumerate(pareto_frontier(archive(alg)))
        delegate_pos > length(delegate_ixs) && break
        @inbounds if delegate_ixs[delegate_pos] == candi_ix
            alg.population[delegate_pos] = candi
            delegate_pos += 1 # move to the next delegate
        end
    end
    # inject mutated archive members
    populate_by_mutants(alg, ndelegates)
    alg.select.size = min(new_popsize, max(3, ceil(Int, alg.τ * new_popsize)))
    alg.last_restart = alg.n_steps
    alg.n_restarts+=1
    return alg
end
