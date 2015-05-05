using Distributions

type BorgMOEA{I <: Individual} <: EvolutionaryAlg
  operators::Vector{GeneticOperator}
  population::Population{I}
  archive::Archive{I}
end

# Take one step of Borg MOEA.
function step(b::BorgMOEA, numInParallel = 32)
  # Update operator probabilities based on archive tag counts
  counts = tagcounts(archive) + 1
  probabilities = counts ./ sum(counts)

  # Select the operators to apply based on their probabilities
  distr = Multinomial(length(probabilities), probabilities)
  ops = rand(distr, numInParallel)
  parents = map(1:numInParallel) do

  end
end