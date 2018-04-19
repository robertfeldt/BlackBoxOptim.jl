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
                recombinate,
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
    # Select the operators to apply based on their probabilities
    recomb_op_ix = rand(alg.recombinate_distr)
    recombinate!(alg, recomb_op_ix, alg.recombinate[recomb_op_ix])
    return alg
end

# "kernel function" of step() that would specialize to given xover operator type
function recombinate!(alg::BorgMOEA, recomb_op_ix::Int, recomb_op::CrossoverOperator)
    # select parents for recombination
    n_parents = numparents(recomb_op)
    # parent indices
    parent_indices = select(alg.select, alg.population,
                            isempty(archive(alg)) ? n_parents : n_parents-1)
    if !isempty(archive(alg))
        # get one parent from the archive and copy it to the fitrst transient member
        arch_ix = transient_range(alg.population)[1]
        alg.population[arch_ix] = archive(alg).frontier[rand_frontier_index(archive(alg))]
        push!(parent_indices, arch_ix)
    end
    # Crossover parents and target
    children = acquire_candis(alg.population, numchildren(recomb_op))
    apply!(recomb_op, Individual[child.params for child in children],
           zeros(Int, length(children)), alg.population, parent_indices)
    for child in children
        child.extra = recomb_op
        child.tag = recomb_op_ix
        process_candidate!(alg, child, parent_indices[1])
    end
end

function process_candidate!(alg::BorgMOEA, candi::Candidate, ref_index::Int)
    apply!(alg.embed, candi.params, alg.population, ref_index)
    reset_fitness!(candi, alg.population)
    ifitness = fitness(update_fitness!(alg.evaluator, candi)) # implicitly updates the archive
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
    alg
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

function update_population_fitness!(alg::BorgMOEA)
    fs = fitness_scheme(alg.evaluator.archive)
    for i in 1:popsize(alg.population)
        if isnafitness(fitness(alg.population, i), fs)
            candi = acquire_candi(alg.population, i)
            update_fitness!(alg.evaluator, candi)
            alg.population.fitness[candi.index] = candi.fitness
            release_candi(alg.population, candi)
        end
    end
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
    mutant
end

function populate_by_mutants(alg::BorgMOEA, last_nonmutant::Int)
    popsz = popsize(alg.population)
    for i in (last_nonmutant+1):popsz
        mutant = acquire_mutant(alg, i, last_nonmutant)
        update_fitness!(alg.evaluator, mutant)
        accept_candi!(alg.population, mutant)
    end
end

"""
Restart Borg MOEA.

Resize and refills the population from the archive.
"""
function restart!(alg::BorgMOEA)
    notify!(archive(alg), :restart)
    frontier_ixs = occupied_frontier_indices(archive(alg))
    narchived = length(frontier_ixs)
    new_popsize = max(alg.min_popsize, ceil(Int, alg.γ * narchived))
    # fill populations with the solutions from the archive
    resize!(alg.population, new_popsize)
    last_archived = min(narchived, new_popsize)
    @inbounds for i in 1:last_archived
        alg.population[i] = archive(alg).frontier[frontier_ixs[i]]
    end
    # inject mutated archive members
    populate_by_mutants(alg, last_archived)
    alg.select.size = min(new_popsize, max(3, ceil(Int, alg.τ * new_popsize)))
    alg.last_restart = alg.n_steps
    alg.n_restarts+=1
    return alg
end
