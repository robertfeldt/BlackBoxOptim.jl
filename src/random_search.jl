type RandomSearcher <: Optimizer
  name::ASCIIString
  search_space::SearchSpace
  best
  best_fitness
  RandomSearcher(searchSpace) = new("Random Search", searchSpace, nothing, nothing)
end

function ask(rs::RandomSearcher)
  # Just randomly generate a new individual and return it with a dummy index
  # (since we do not have a population there are no indices).
  [(rand_individual(rs.search_space), 1)]
end

function tell!(rs::RandomSearcher, rankedCandidates)
  candidate, index, fitness = rankedCandidates[1]
  if rs.best == nothing || fitness < rs.best_fitness
    rs.best = candidate
    rs.best_fitness = fitness
    return 1
  else
    return 0
  end
end

function random_search(parameters)
  RandomSearcher(parameters[:SearchSpace])
end