facts("Parallel population optimizer") do

@everywhere using BlackBoxOptim

context("optimizing small problem") do
  rosenbrock2d(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

  res = bboptimize(rosenbrock2d; Method = :parallel_population_optimizer,
      SearchSpace = [(-5.0, 5.0), (-2.0, 2.0)], MaxTime = 100.0,
      ShowTrace = true, MigrationSize = 2, MigrationPeriod = 100)
  @fact size(best_candidate(res)) => (2,)
  @fact typeof(best_fitness(res)) => Float64
  @fact best_fitness(res) => less_than(100.0)
  pop = population(res)
  @fact numdims(pop) --> 2
  @fact popsize(pop) > 0 --> true
end

context("propagating exceptions from workers to the master") do
  # 0.01 chance to get domain error
  bogus(x) = sqrt(x[1]-0.01)
  # 0.01 chance to raise exception prematurely during bbsetup()
  res = bbsetup(bogus; Method = :parallel_population_optimizer,
      SearchSpace = [(0.0, 1.0)], MaxTime = 100.0,
        ShowTrace = true, MigrationSize = 2, MigrationPeriod = 100)
  @fact_throws RemoteException bboptimize(res)
end

#= doesn't work with PopulationMatrix
context("worker method that uses PopulationMatrix") do
  rosenbrock2d(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

  res = bboptimize(rosenbrock2d; Method = :parallel_population_optimizer,
      WorkerMethod = :separable_nes,
      SearchSpace = [(-5.0, 5.0), (-2.0, 2.0)], MaxTime = 100.0,
      ShowTrace = true, MigrationSize = 2, MigrationPeriod = 100)
  @fact size(best_candidate(res)) => (2,)
  @fact typeof(best_fitness(res)) => Float64
  @fact best_fitness(res) => less_than(100.0)
  pop = population(res)
  @fact numdims(pop) --> 2
  @fact popsize(pop) > 0 --> true
end
=#
end
