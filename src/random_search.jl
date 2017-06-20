"""
  Optimize by randomly generating the candidates.
"""
type RandomSearcher{S<:SearchSpace} <: AskTellOptimizer
  name::String
  search_space::S
  best_fitness          # FIXME fitness type should be known
  best::Individual

  RandomSearcher{S}(searchSpace::S) where S = new("Random Search", searchSpace, nothing)
end

RandomSearcher{S<:SearchSpace}(searchSpace::S) = RandomSearcher{S}(searchSpace)

function ask(rs::RandomSearcher)
  # Just randomly generate a new individual and return it with a dummy index
  # (since we do not have a population there are no indices).
  Candidate{Float64}[Candidate{Float64}(rand_individual(rs.search_space), 1)]
end

function tell!{F}(rs::RandomSearcher, rankedCandidates::Vector{Candidate{F}})
  candidate = rankedCandidates[1]
  if rs.best_fitness == nothing || candidate.fitness < rs.best_fitness
    rs.best = candidate.params
    rs.best_fitness = candidate.fitness
    return 1
  else
    return 0
  end
end

function random_search(problem::OptimizationProblem, parameters::Parameters)
  RandomSearcher(search_space(problem))
end

function random_search(ss::SearchSpace)
  RandomSearcher(ss)
end
