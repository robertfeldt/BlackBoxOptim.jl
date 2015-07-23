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
    println("New RW iteration")
    step!(rw)
    sleep(0.5 * rand())
  end
end

history = Input( (Int, Float64)[] )

function report_new_min(i::Int, newmin::Float64)
  push!(history, push!(value(history), (i, newmin)))
end

rw = RandWalker(report_new_min)
t = @async rand_walk!(rw)

using Gadfly

function main(window)
    push!(window.assets, "widgets")

    lift(history) do  hist
      xvals = map(first, hist)
      yvals = map(t -> t[2], hist)

      println("Changed!, hist = ", hist)

        vbox(
          h1("Fitness progress") |> emph,
          hbox("Best fitness: ", hskip(1em), @sprintf("%.3e", yvals[end]) |> emph),
          hbox("Found at iteration: ", hskip(1em), string(xvals[end]) |> emph),
          # Scale.x_log10, Stat.step
          plot(x = xvals, y = yvals, Geom.point, Geom.line) |> drawing(8inch, 6inch)
        ) |> packacross(center)
    end
end
