using HTTP, JSON

const VegaLiteWebsocketFrontEndTemplate = """
<!DOCTYPE html>
<head>
  <meta charset="utf-8">
  <script src="https://cdn.jsdelivr.net/npm/vega@3"></script>
	<script src="https://cdn.jsdelivr.net/npm/vega-lite@2"></script>
	<script src="https://cdn.jsdelivr.net/npm/vega-embed@3"></script>
</head>

<body>
  <div id="vis"></div>

  <script>
    const vegalitespec = %%VEGALITESPEC%%;
    vegaEmbed('#vis', vegalitespec, {defaultStyle: true})
    .then(function(result) {
        const view = result.view;
        const port = %%SOCKETPORT%%;
        const conn = new WebSocket("ws://127.0.0.1:" + port);

        conn.onopen = function(event) {
          // insert data as it arrives from the socket
          conn.onmessage = function(event) {
            console.log(event.data);
            // Use the Vega view api to insert data
            var newentries = JSON.parse(event.data);
            view.insert("table", newentries).run();
          }
        }
      })
    .catch(console.warn);
  </script>
</body>
"""

const VegaLiteFitnessOverTimeSpecTemplate = """
{
    "\$schema": "https://vega.github.io/schema/vega-lite/v4.json",
    "description": "Fitness value over time",
    "width": %%WIDTH%%,
    "height": %%HEIGHT%%,
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
            "field": "Fitness",
            "type": "quantitative",
            "scale": {"type": "log"}
        }
    }
}
"""

rand_websocket_port() = 9000 + rand(0:42)

function vegalite_websocket_frontend(vegalitespec::String = VegaLiteFitnessOverTimeSpecTemplate;
    width::Int = 800, height::Int = 600, port::Int = rand_websocket_port())
    res = replace(VegaLiteWebsocketFrontEndTemplate, "%%VEGALITESPEC%%" => vegalitespec)
    res = replace(res, "%%WIDTH%%" => string(width))
    res = replace(res, "%%HEIGHT%%" => string(height))
    replace(res, "%%SOCKETPORT%%" => string(port))
end

function serve_html_file(file::String)
    HTTP.serve() do request::HTTP.Request
        try
            return HTTP.Response(file)
        catch e
            return HTTP.Response(404, "Error: $e")
        end
    end
end

# Serve a VegaLite front-end on html and then push updates to data
# over a websocket so that the VegaLite graph is updated.
mutable struct VegaLiteGraphFitnessGraph
    verbose::Bool
    httpport::Int
    websocketport::Int
    mindelay::Float64
    data::Vector{Any}
    last_sent_index::Int
    starttime::Float64
    function VegaLiteGraphFitnessGraph(verbose::Bool = true,
        httpport::Int = 8081, websocketport::Int = 9000, mindelay::Float64 = 1.0)
        @assert mindelay > 0.0
        new(verbose, httpport, websocketport, mindelay, Any[], 0, 0.0)
    end
end

timestamp(t = time()) = Libc.strftime("%Y-%m-%d %H:%M.%S", t)
log(vlg::VegaLiteGraphFitnessGraph, msg) = vlg.verbose ? println(timestamp(), ": ", msg) : nothing

function add_data!(vlg::VegaLiteGraphFitnessGraph, newentry::Dict)
    if length(vlg.data) < 1
        vlg.starttime = time()
    end
    if !haskey(newentry, "Time")
        newentry["Time"] = time() - vlg.starttime
    end
    log(vlg, "Adding data $newentry")
    push!(vlg.data, newentry)
end

function serve(vlg::VegaLiteGraphFitnessGraph)
    @async websocket_push_data_loop(vlg.websocketport, vlg.mindelay) do ws
        send_new_data_on_socket(vlg, ws)
    end
    frontend_html = vegalite_websocket_frontend(; port = vlg.websocketport)
    @async serve_html_file(frontend_html)
    println("Serving frontend on http://127.0.0.1:$(vlg.httpport)")
    sleep(3.0) # To give people time to copy the url...
end

hasnewdata(vlg::VegaLiteGraphFitnessGraph) = length(vlg.data) > vlg.last_sent_index

function send_new_data_on_socket(vlg::VegaLiteGraphFitnessGraph, ws)
    if hasnewdata(vlg)
        len = length(vlg.data)
        newdata = vlg.data[(vlg.last_sent_index+1):len]
        log(vlg, "Sending data $newdata")
        write(ws, JSON.json(newdata))
        vlg.last_sent_index = len
    end
end

function websocket_push_data_loop(fn::Function, wsport::Int, mindelay = 1.0)
    HTTP.WebSockets.listen("127.0.0.1", UInt16(wsport)) do ws
        while true
            fn(ws) # Call with websocket
            sleep(mindelay + rand())
        end
    end
end
