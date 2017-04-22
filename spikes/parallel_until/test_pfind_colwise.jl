start = time()
N = 100
Rows = 3
MaxValue = 10
#MaxNumProcs = 4
Reps = 10

@everywhere delay = 0.01

#times = zeros(Reps, MaxNumProcs)
times = zeros(Reps, 1)

#for numproc in 1:MaxNumProcs
  numproc = nprocs()
  #println("With $(numproc) procs")

  require("pfind_colwise.jl") # Reload to ensure it is on all procs
 
  @everywhere mycond(column) = fake_slow_condition(column, delay, 30)

  for rep in 1:Reps
    #println("  Rep $rep")
    array = ifloor(MaxValue * rand(Rows, N))

    # We will be looking for the first column that sums to 30. Make sure it exists.
    suma = sum(array, 1)
    i = findfirst(suma, 30)
    if i == 0
      # Put it somewhere around the middle
      #i = ifloor(rand(0.45N:0.55*N))
      i = rand(1:N)
      array[:,i] = MaxValue * ones(Int, Rows, 1)
    end

    tic()
      idx, res = pfind_colwise(mycond, array)
    #times[rep, numproc] = toq()
    times[rep] = toq()
    @assert idx == i
  end

  #addprocs(1) # Add one proc for next round
#end

#println(mean(times[2:Reps,:], 2))
mean_time = mean(times[2:Reps,:])
println("Nprocs = $(nprocs()), Delay = $delay, Avg time = $(mean_time), Total time = $(time() - start)")