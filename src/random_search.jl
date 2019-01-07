"""
Optimize by randomly generating the candidates.
"""
mutable struct RandomSearcher{SS<:SearchSpace} <: AskTellOptimizer
    name::String
    search_space::SS
    best_fitness          # FIXME fitness type should be known
    best::Individual

    RandomSearcher(search_space::SS) where {SS<:SearchSpace} =
        new{SS}("Random Search", search_space, nothing)
end

function ask(rs::RandomSearcher)
    # Just randomly generate a new individual and return it with a dummy index
    # (since we do not have a population there are no indices).
    [Candidate{Float64}(rand_individual(rs.search_space), 1)]
end

function tell!(rs::RandomSearcher, rankedCandidates::Vector{Candidate{F}}) where F
    candidate = rankedCandidates[1]
    if rs.best_fitness == nothing || candidate.fitness < rs.best_fitness
        rs.best = candidate.params
        rs.best_fitness = candidate.fitness
        return 1
    else
        return 0
    end
end

random_search(problem::OptimizationProblem, parameters::Parameters) =
    RandomSearcher(search_space(problem))

random_search(search_space::SearchSpace) =
    RandomSearcher(search_space)
