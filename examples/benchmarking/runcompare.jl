scripts_dir = @__DIR__
reports_dir = (length(ARGS) > 0 ? ARGS[1] : "examples/benchmarking")

isdir(reports_dir) || throw(error("Benchmark results folder $reports_dir not found"))

bfiles = filter!(f->occursin(r"^benchmark_runs_\d.+\.csv$", f), readdir(reports_dir))

isempty(bfiles) && throw(error("No benchmark runs found"))
bf = joinpath(reports_dir, bfiles[end])

println("Comparing to benchmark: ", bf)

nreps = (length(ARGS) > 1 ? parse(Int, ARGS[2]) : 10)


@info @__DIR__

cmd = `$(Sys.BINDIR)/julia --color=yes -L $scripts_dir/activate_pwd_env.jl -- $scripts_dir/compare_optimizers.jl compare --benchmarkfile $bf --numreps $nreps`

run(cmd)
