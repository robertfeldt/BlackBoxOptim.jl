# Idea is that instead of having our own Tracer functionality we would use the
# new Logging interface in Julia StdLib from 0.7. Here we test how this could work.
using Logging: AbstractLogger, with_logger, global_logger

timestamp(t = time()) = Libc.strftime("%Y-%m-%d %H:%M.%S", t)
elapsed(st, t = time()) = round(t - st; digits=4)

# Try the default logger
@info "Optimization step" time=timestamp() elapsed=elapsed(ST) best_fitness=3.4

# Let's make a CSVLogger to log to CSV file.
struct CSVLogger <: AbstractLogger
    starttime::Float64
    filename::String
    filehandle::IOStream
    columns::Vector{Symbol}
    group::Union{Nothing, AbstractString, Vector{AbstractString}}
end
function CSVLogger(filename::String, columns::Vector{Symbol}, 
    group::Union{Nothing, AbstractString, Vector{AbstractString}} = nothing)
    if isfile(filename)
        error("Cannot open csv log file $filename, since it already exists")
    end
    fh = open(filename, "w")
    writeheader(fh, columns)
    CSVLogger(time(), filename, fh, columns, group)
end

function writeheader(fh::IOStream, columnnames::Vector{Symbol})
    names = vcat(["Time", "Elapsed"], map(string, columnnames))
    println(fh, join(names, ","))
    flush(fh)
end

reset!(l::CSVLogger) = (l.starttime = time())

import Logging: min_enabled_level, shouldlog, handle_message
@inline min_enabled_level(l::CSVLogger) = min_enabled_level(global_logger())

function shouldlog(l::CSVLogger, level, _module, group, id)
    @show (group, typeof(group))
    true
end

@inline shouldlog_message(g::Nothing, message::AbstractString) = true # Log everything if no group set in CSVLogger.group
@inline shouldlog_message(g::AbstractString, message::AbstractString) = (g == message)
@inline shouldlog_message(g::Vector{AbstractString}, message::AbstractString) = in(message, g)

function handle_message(l::CSVLogger, level, message, _module, group, id, file, line; kwargs...)
    t = time()
    ts = timestamp(t)
    el = elapsed(l.starttime, t)
    if shouldlog_message(l.group, message)
        @show (ts, el, level, message, _module, group, id, file, line, kwargs)
        values = map(k -> kwargs[k], l.columns)
        @show values
        println(l.filehandle, ts * "," * string(el) * "," * join(values, ","))
        flush(l.filehandle)
    end
end

l = CSVLogger("t.csv", [:BestFitness]) # We will only log best_fitness (and add time and elapsed num seconds)
with_logger(l) do
    @info "Optimization step" BestFitness=3.4 Fitness=45.6 Ind="1+3"
end

with_logger(l) do
    @info "Optimization step" BestFitness=2.9714 Fitness=32.6 Ind="1+578"
end