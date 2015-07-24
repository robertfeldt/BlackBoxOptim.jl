type RandWalker
  i::Int
  x::Float64
  y::Float64
  ymin::Float64
  reporter
  function RandWalker(reporter = nothing)
    new(0, 0.0, 0.0, 0.0, reporter)
  end
end

function step!(rw::RandWalker)
  rw.i += 1
  rw.x += randn()
  rw.y += randn()
  if rw.y < rw.ymin
    println("New min found: ($(rw.i),$(rw.y))")
    rw.ymin = rw.y
    if rw.reporter != nothing
      rw.reporter(rw.i, rw.ymin)
    end
  end
end

function rand_walk!(rw::RandWalker)
  while true
    print("."); flush(STDOUT);
    step!(rw)
    sleep(1.0 * rand())
  end
end

fitness_history = Input( (Int, Float64)[] )

function report_new_min(i::Int, newmin::Float64)
  push!(fitness_history, push!(value(fitness_history), (i, newmin)))
end

rw = RandWalker(report_new_min)
t = @async rand_walk!(rw)

using Gadfly

function fitness_plot(iterations, fitnesses)
  plot(x = iterations, y = fitnesses,
    Geom.point, Geom.line,
    Guide.xlabel("Iteration"), Guide.ylabel("Fitness")
  ) |> drawing(10inch, 6inch)
end

function main(window)
    push!(window.assets, "widgets")
    local content

    lift(fitness_history) do hist
      println("Changed! len = ", length(hist))
      if length(hist) >= 1
        its  = map(first, hist)
        fits = map(t -> t[2], hist)
        content = vbox(
          hbox("Best fitness:       ", hskip(1em), @sprintf("%.3e", fits[end]) |> emph),
          hbox("Found at iteration: ", hskip(1em), string(its[end]) |> emph),
          fitness_plot(its, fits),
        )
      else
        content = hbox("No min found yet...")
      end

      vbox(
        h1("Fitness progress") |> emph,
        content
      ) |> packacross(center)
    end
end
