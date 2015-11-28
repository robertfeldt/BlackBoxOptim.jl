using BlackBoxOptim
using FactCheck
using Compat

NumTestRepetitions = 100

if nprocs() < 4
  addprocs(4-nprocs()) # add procs for parallel population optimizer
end
@everywhere using BlackBoxOptim
