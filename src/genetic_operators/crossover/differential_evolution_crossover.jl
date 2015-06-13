abstract DiffEvoCrossoverOperator{NP,NC} <: CrossoverOperator{NP,NC}

type DiffEvoRandBin1{T <: Real} <: DiffEvoCrossoverOperator{3,1}
  f::Float64
  cr::Float64
  DiffEvoRandBin1(f = 0.65, cr = 0.4) = new(f, cr)
end

function apply{T <: Real}(xo::DiffEvoRandBin1{T}, p1::Vector{T}, p2::Vector{T}, p3::Vector{T})
  donor = p3 .+ xo.f .* (p1 .- p2)
  de_crossover_binomial(p3, donor)
end

type DiffEvoRandBin2{T <: Real} <: DiffEvoCrossoverOperator{5,1}
  f::Float64
  cr::Float64
  DiffEvoRandBin2(f = 0.65, cr = 0.4) = new(f, cr)
end

function apply{T <: Real}(xo::DiffEvoRandBin2{T}, p1::Vector{T}, p2::Vector{T}, p3::Vector{T}, p4::Vector{T}, p5::Vector{T})
  donor = p3 .+ de.f * (p1 .- p2) .+ de.f * (p4 .- p5)
  de_crossover_binomial(p3, donor)
end

function de_crossover_binomial{T <: Real}(target::Vector{T}, donor::Vector{T})
  trial = copy(target)

  # Always ensure at least one value from donor is copied to trial vector.
  jrand = rand(1:length(trial))
  trial[jrand] = donor[jrand]

  # Now crossover randomly for the rest of the indices.
  switch = rand(length(trial)) .<= de.cr
  trial[switch] = donor[switch]

  return trial
end
