# simple random selector
type SimpleSelector <: IndividualsSelector
end

function select(::SimpleSelector, population, numSamples::Int)
  sample(1:popsize(population), numSamples; replace = false)
end
