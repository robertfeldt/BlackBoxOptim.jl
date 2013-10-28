function rosenbrock2d(x)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

facts("bboptimize") do
  best, fitness = bboptimize(rosenbrock2d, (-5.0, 5.0); dimensions = 2)
  @fact fitness < 0.01 => true
end