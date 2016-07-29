using BlackBoxOptim

if length(ARGS) >= 1
  maxtime = parse(Float64, ARGS[1])
else
  maxtime = 0.3 # Time used in orignal (flaky) tests...
end

if length(ARGS) >= 2
  N = parse(Int, ARGS[2])
else
  N = 30
end

schaffer1(x) = (sumabs2(x), sumabs2(x .- 2.0))

results1 = Float64[]
results2 = Float64[]

totaltime = 0.0

for i in 1:(N+1)
  ctrl = bbsetup(schaffer1; Method = :borg_moea,
    SearchRange = [(-10.0, 10.0), (-10.0, 10.0)], TraceMode = :silent,
    FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true), ϵ=0.01)
  tic()
  @time res = bboptimize(ctrl, MaxTime = maxtime)
  t = toq()
  if i > 1 # Skip first since there is compile time involved
    totaltime += t
    f = best_fitness(res)
    push!(results1, f[1])
    push!(results2, f[2])
  end
end

function report(vs, label, t = nothing)
  r(x) = round(x, 3)
  println(label, ": ", r(mean(vs)), " +/- ", r(std(vs)), 
    " (", r(minimum(vs)), "-", r(maximum(vs)), ")")
  if t != nothing
    println("Time: ", r(t))
  end
end

report(results1, "Fitness 1")
report(results2, "Fitness 2", totaltime)

# Results below from running on Robert's MacBook Pro 13" from 2015, on 2016-07-29:

# julia04 -L src/BlackBoxOptim.jl examples/stability_test.jl 0.3 50
# Fitness 1: 1.999 +/- 0.007 (1.974-2.016)
# Fitness 2: 2.001 +/- 0.007 (1.985-2.026)
# Time: 15.028

# julia05 -L src/BlackBoxOptim.jl examples/stability_test.jl 0.3 50
# Fitness 1: 2.031 +/- 0.234 (1.655-2.85)
# Fitness 2: 2.059 +/- 0.423 (1.48-4.653)
# Time: 15.112

# julia04 -L src/BlackBoxOptim.jl examples/stability_test.jl 3.0 50
# Fitness 1: 2.0 +/- 0.005 (1.995-2.005)
# Fitness 2: 2.0 +/- 0.005 (1.995-2.005)
# Time: 147.778

# julia05 -L src/BlackBoxOptim.jl examples/stability_test.jl 3.0 50
# Fitness 1: 2.004 +/- 0.019 (1.955-2.048)
# Fitness 2: 1.996 +/- 0.019 (1.952-2.046)
# Time: 151.186