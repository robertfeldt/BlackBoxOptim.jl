# SBX = Simulated Binary Crossover

calc_eta_exponent(eta) = 1 / (eta + 1)

type SimulatedBinaryCrossover <: CrossoverOperator
  eta
  etaexponent # Precalc the eta exponent used in calc of beta
  SimulatedBinaryCrossover(eta = 3.0) = new(eta, calc_eta_exponent(eta))
end

function apply{T <: Real}(sbx::SimulatedBinaryCrossover, p1::T, p2::T)
  u = rand()
  if u < 0.5
    beta = (2*u)^sbx.etaexponent
  else
    beta = (1/(2*(1-u)))^sbx.etaexponent
  end
  c1 = 0.5 * ( (1+beta) * p1 + (1-beta) * p2 )
  c2 = 0.5 * ( (1-beta) * p1 + (1+beta) * p2 )
  return (c1, c2)
end
