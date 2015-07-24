type RandWalker
  i::Int
  x::Float64
  y::Float64
  ymin::Float64
  reporter
  function RandWalker(reporter = nothing)
    new(0, 0.0, 0.0, Inf, reporter)
  end
end

function step!(rw::RandWalker)
  rw.i += 1
  rw.x += randn()
  rw.y += (randn() - 0.05) # Slight downward trend to make more fitness improvement happen...
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
    sleep(0.4 * rand())
  end
end

using Gadfly

function fitness_plot(iterations, fitnesses)
  plot(x = iterations, y = fitnesses,
    Geom.point, Geom.line,
    Guide.xlabel("Iteration"), Guide.ylabel("Fitness"),
    Stat.step,
    #Scale.y_log10,
  ) |> drawing(10inch, 6inch)
end

function main(window)
  push!(window.assets, "widgets")

  # If we don't already have a fitness history we create one and start a
  # random walker that updates it.
  if !haskey(task_local_storage(), :fhistory)
    println("Starting random walker...")
    fitness_history = Input( (Int, Float64)[] )

    function report_new_min(i::Int, newmin::Float64)
      push!(fitness_history, push!(value(fitness_history), (i, newmin)))
    end

    # Create random walker, start task for walking
    rw = RandWalker(report_new_min)
    t = @async rand_walk!(rw)

    task_local_storage(:fhistory, fitness_history)
    task_local_storage(:rwtask, t)
  end

  local content

  lift(task_local_storage(:fhistory)) do hist
    println("Changed! len = ", length(hist))

    if length(hist) >= 1
      its  = map(first, hist)
      fits = map(t -> t[2], hist)
      content = vbox(
        hbox("Best fitness:       ", hskip(1em), @sprintf("%.3f", fits[end]) |> emph),
        hbox("Found at iteration: ", hskip(1em), string(its[end]) |> emph),
        fitness_plot(its, fits),
      ) |> packacross(center)
    else
      content = hbox("No min found yet...")
    end

    task_started = task_running = false
    if haskey(task_local_storage(), :rwtask)
      t = task_local_storage(:rwtask)
      task_running = !istaskdone(t)
      #task_started = istaskstarted(t) # Only in 0.4
    end

    vbox(
      h1("Fitness progress") |> emph,
      content,
      #hbox("Task started: ", string(task_started)), # Only in 0.4
      #hbox("Task running: ", hskip(1em), string(task_running)),
    ) |> packacross(center)
  end
end
