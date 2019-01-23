using HttpServer, WebSockets, JSON

timestamp(t = time()) = Libc.strftime("%Y-%m-%d %H:%M.%S", t)

function start_http_and_websocket_servers(httpcontent, httpport = 8081, wsport = 8083)
    http = HttpHandler() do req::Request, res::Response
        @show req
        Response( servedcontent )
    end
    httpserver = Server( http )
    println("Starting http server on port $httpport")
    @async run(httpserver, httpport)

    wsh = WebSocketHandler() do req,client
        println("WebSocket opened from $client:\n $req")

        # Send meta-info
        starttime = time()
        ts = timestamp(starttime)
        metainfo = Dict(:info => "meta", :StartTime => ts)
        write(client, JSON.json(metainfo))

        # For sending optimization events we do not set the info key. This saves some time
        # in the json conversion on both sides. We reuse one and the same Dict to save on mem.
        obj = Dict(:ElapsedSeconds => 0, :BestFitness => 1000 + rand()*1000, :NumEvals => 1)

        while true
            sleep(0.1 + rand() * 0.3)
            obj[:ElapsedSeconds] = t = time() - starttime
            obj[:BestFitness] -= (rand() * 100 / (3 + log10(100 + 0.1 * t)))
            obj[:NumEvals] += 1
            try
                write(client, JSON.json(obj))
            catch err
                @show err
                break
            end
        end
    end

    wsserver = Server(wsh)
    println("Starting web socket server on port $wsport")
    @async run(wsserver, wsport)

    return httpserver, wsserver
end

servedcontent = read(open("vega_lite_updating_from_websocket.html"), String)
start_http_and_websocket_servers(servedcontent, 8081, 8083)
sleep(120.0)