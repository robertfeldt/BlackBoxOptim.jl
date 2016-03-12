nbfile = Libc.strftime("benchmark_runs_%y%m%d.csv", time())

if isfile(nbfile)
  println("Deleting existing file: ", nbfile)
  rm(nbfile)
end

nreps = (length(ARGS) > 0 ? parse(Int, ARGS[1]) : 10)

cmd = `julia --color=yes -L ../../src/BlackBoxOptim.jl compare_optimizers.jl update --benchmarkfile $nbfile --numreps $nreps`

run(cmd)