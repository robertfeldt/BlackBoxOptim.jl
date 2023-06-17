replace_template_param(template::AbstractString, param_to_value::Pair{Symbol, <:Any}) =
    replace(template, string("%%", first(param_to_value), "%%") => string(last(param_to_value)))

"""
Specification of a plot for the real-time tracking of the fitness progress.

To use the [*VegaLite*](https://vega.github.io/vega-lite/) front-end via
*BlackBoxOptimRealtimePlotServerExt* extension, *HTTP.jl* and *Sockets.jl* are required.
"""
mutable struct RealtimePlot{E}
    spec::String
    verbose::Bool
    data::Vector{Any}
    last_sent_index::Int
    starttime::Float64
    stoptime::Union{Float64, Nothing}

    function RealtimePlot{E}(template::AbstractString;
        verbose::Bool = false,
        spec_kwargs...
    ) where E
        @assert E isa Symbol
        spec = reduce(replace_template_param, spec_kwargs,
                      init = template)
        new{E}(spec, verbose, Any[], 0, 0.0, nothing)
    end
end

timestamp(t = time()) = Libc.strftime("%Y-%m-%d %H:%M.%S", t)
printmsg(plot::RealtimePlot, msg) = plot.verbose ? println(timestamp(), ": ", msg) : nothing

function shutdown!(plot::RealtimePlot)
    plot.stoptime = time()
end

function Base.push!(plot::RealtimePlot, newentry::AbstractDict)
    if length(plot.data) < 1
        plot.starttime = time()
    end
    if !haskey(newentry, "Time")
        newentry["Time"] = time() - plot.starttime
    end
    printmsg(plot, "Adding data $newentry")
    push!(plot.data, newentry)
end

hasnewdata(plot::RealtimePlot) = length(plot.data) > plot.last_sent_index

const VegaLiteMetricOverTimePlotTemplate = """
{
    "\$schema": "https://vega.github.io/schema/vega-lite/v4.json",
    "description": "%%metric%% value over time",
    "width": %%width%%,
    "height": %%height%%,
    "padding": {"left": 20, "top": 10, "right": 10, "bottom": 20},
    "data": {
      "name": "table"
    },
    "mark": "line",
    "encoding": {
        "x": {
            "field": "Time",
            "type": "quantitative"
        },
        "y": {
            "field": "%%metric%%",
            "type": "quantitative",
            "scale": {"type": "log"}
        }
    }
}
"""

VegaLiteMetricOverTimePlot(; metric::String = "Fitness",
                           width::Integer = 800, height::Integer = 600,
                           kwargs...) =
    RealtimePlot{:VegaLite}(VegaLiteMetricOverTimePlotTemplate; metric, width, height, kwargs...)

"""
    fitness_plot_callback(plot::RealtimePlot, oc::OptRunController)

[OptController](@ref) callback function that updates the real-time fitness plot.
"""
function fitness_plot_callback(plot::RealtimePlot, oc::OptRunController)
    push!(plot, Dict("num_steps" => num_steps(oc),
                     "Fitness" => best_fitness(oc)))
    if oc.stop_reason != ""
        @info "Shutting down realtime plot"
        shutdown!(plot)
    end
end
