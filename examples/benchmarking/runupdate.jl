scripts_dir = @__DIR__
reports_dir = (length(ARGS) > 0 ? ARGS[1] : "examples/benchmarking")

nbfile = joinpath(reports_dir, Libc.strftime("benchmark_runs_%y%m%d.csv", time()))

if isfile(nbfile)
    @info("Deleting existing benchmark file: ", nbfile)
    rm(nbfile)
end

nreps = (length(ARGS) > 0 ? parse(Int, ARGS[2]) : 10)

cmd = `$(Sys.BINDIR)/julia --color=yes -L $scripts_dir/activate_pwd_env.jl -- $scripts_dir/compare_optimizers.jl update --benchmarkfile $nbfile --numreps $nreps`

run(cmd)
