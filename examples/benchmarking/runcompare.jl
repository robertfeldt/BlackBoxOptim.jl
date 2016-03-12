bfiles = filter(f->typeof(match(r"^benchmark_runs_\d.+\.csv$", f))!=Void, readdir("."))

bf = bfiles[end]

println("Comparing to benchmark: ", bf)

nreps = (length(ARGS) > 0 ? parse(Int, ARGS[1]) : 10)

cmd = `julia --color=yes -L ../../src/BlackBoxOptim.jl compare_optimizers.jl compare --benchmarkfile $bf --numreps $nreps`

run(cmd)