abstract OptimizationResults{F}

immutable SimpleOptimizationResults{F,RT} <: OptimizationResults{F}
  method::ASCIIString # Symbol instead or flexible?
  best_fitness::F
  best_candidate::RT
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
# FIXME should be it be enabled only for MinimizingFitnessScheme?
Base.minimum(or::OptimizationResults) = best_candidate(or)
f_minimum(or::OptimizationResults) = best_fitness(or)
# FIXME lookup stop_reason
iteration_converged(or::OptimizationResults) = iterations(or) >= parameters(or)[:MaxSteps]
