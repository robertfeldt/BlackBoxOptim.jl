# Polynomial mutation as presented in the paper:
#  Deb and Deb (2012), "Analyzing Mutation Schemes for Real-Parameter Genetic Algorithms"

abstract GeneticOperator
abstract MutationOperator <: GeneticOperator

type PolynomialMutation <: MutationOperator
  lowbound
  highbound
  eta
  PolynomialMutation(lowbound = 0.0, highbound = 1.0, eta = 100.0) = new(lowbound, highbound, eta)
end

# Apply polynomial mutation to one value.
function apply(pm::PolynomialMutation, value::Float64)
  u = rand()
  if u <= 0.5
    deltaL = (2*u)^(1/(1+pm.eta)) - 1
    return (value + deltaL * (value - pm.lowbound))
  else
    deltaR = 1 - (2*(1-u))^(1/(1+pm.eta))
    return (value + deltaR * (pm.highbound - value))
  end
end