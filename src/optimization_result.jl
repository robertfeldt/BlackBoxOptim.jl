abstract OptimizationResults

immutable SingleObjectiveOptimizationResults{RT,FT} <: OptimizationResults
  method::ASCIIString # Symbol instead or flexible?
  best_candidate::RT
  best_fitness::FT
  stop_reason::ASCIIString # FIXME turn into type hierarchy of immutable reasons with their attached info
  iterations::Int
  start_time::Float64
  elasped_time::Float64 # seconds
  parameters::Parameters
  f_calls::Int
end

best_candidate(or::OptimizationResults) = or.best_candidate
best_fitness(or::OptimizationResults) = or.best_fitness
stop_reason(or::OptimizationResults) = or.stop_reason
iterations(or::OptimizationResults) = or.iterations
start_time(or::OptimizationResults) = or.start_time
elapsed_time(or::OptimizationResults) = or.elapsed_time
parameters(or::OptimizationResults) = or.parameters
f_calls(or::OptimizationResults) = or.f_calls

# Alternative nomenclature that mimics Optim.jl more closely.
import Base.minimum
Base.minimum(or::OptimizationResults) = or.best_candidate
f_minimum(or::OptimizationResults) = or.best_fitness
iteration_converged(or::OptimizationResults) = iterations(or) >= parameters(or)[:MaxSteps]
