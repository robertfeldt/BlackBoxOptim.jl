abstract DiffEvoCrossoverOperator{NP,NC} <: CrossoverOperator{NP,NC}

# FIXME is it possible somehow to do arithmetic operations with N?
immutable DiffEvoRandBin{N} <: DiffEvoCrossoverOperator{N,1}
end

typealias DiffEvoRandBin1 DiffEvoRandBin{3}
typealias DiffEvoRandBin2 DiffEvoRandBin{5}

function apply!(xo::DiffEvoRandBin{3}, cr::Real, f::Number, target, pop, parentIndices)
  @assert length(parentIndices) == 3
  p1ix, p2ix, p3ix = parentIndices
  # Always ensure at least one parameter is xovered
  mut_ix = rand(1:length(target))
  for i in 1:length(target)
    if i == mut_ix || rand() <= cr
      target[i] = pop[i,p3ix] + f .* (pop[i,p1ix] .- pop[i,p2ix])
    end
  end
  return target
end

function apply!(xo::DiffEvoRandBin{5}, cr::Real, f::Number, target, pop, parentIndices)
  @assert length(parentIndices) == 5
  p1ix, p2ix, p3ix, p4ix, p5ix = parentIndices
  # Always ensure at least one parameter is xovered
  mut_ix = rand(1:length(target))
  for i in 1:length(target)
    if i == mut_ix || rand() <= cr
      target[i] = pop[i,p3ix] +
                f .* (pop[i,p1ix] .- pop[i,p2ix]) +
                f .* (pop[i,p4ix] .- pop[i,p5ix])
    end
  end
  return target
end
