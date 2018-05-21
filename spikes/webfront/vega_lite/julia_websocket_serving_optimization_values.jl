using HttpServer
using WebSockets
using JSON

timestamp(t = time()) = Libc.strftime("%Y-%m-%d %H:%M.%S", t)

wsh = WebSocketHandler() do req,client
    println("WebSocket opened from $client:\n $req")

    # Send meta-info
    starttime = time()
    ts = timestamp(starttime)
    @show ts
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

server = Server(wsh)
run(server,8082)