module BlackBoxOptimRealtimePlotServerExt

using HTTP, Sockets, JSON
using BlackBoxOptim: RealtimePlot, replace_template_param, hasnewdata, printmsg

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

rand_websocket_port() = 9000 + rand(0:42)

frontend_html(vegalitespec::String, socketport::Integer) =
    reduce(replace_template_param, [
        :SOCKETPORT => string(socketport),
        :VEGALITESPEC => vegalitespec,
    ], init = VegaLiteWebsocketFrontEndTemplate)

function static_content_handler(content::AbstractString, request::HTTP.Request)
    try
        return HTTP.Response(content)
    catch e
        return HTTP.Response(404, "Error: $e")
    end
end

function HTTP.serve(plot::RealtimePlot{:VegaLite},
    host=Sockets.localhost, port::Integer = 8081;
    websocketport::Integer = rand_websocket_port(),
    mindelay::Number = 1.0,
    kwargs...
)
    @assert mindelay > 0.0
    @async websocket_serve(plot, websocketport, mindelay)
    printmsg(plot, "Serving VegaLite frontend on http://$(host):$(port)")
    return HTTP.serve(Base.Fix1(static_content_handler, frontend_html(plot.spec, websocketport)),
                      host, port; kwargs...)
end

HTTP.serve!(plot::RealtimePlot, args...; kwargs...) =
    @async(HTTP.serve(plot, args...; kwargs...))

function websocket_serve(plot::RealtimePlot, port::Integer, mindelay::Number)
    HTTP.WebSockets.listen(Sockets.localhost, UInt16(port)) do ws
        while true
            if hasnewdata(plot)
                len = length(plot.data)
                newdata = plot.data[(plot.last_sent_index+1):len]
                printmsg(plot, "Sending data $newdata")
                HTTP.WebSockets.send(ws, JSON.json(newdata))
                plot.last_sent_index = len
            end
            isnothing(plot.stoptime) || break
            sleep(mindelay + rand())
        end
    end
    printmsg(plot, "RealtimePlot websocket stopped")
end

end