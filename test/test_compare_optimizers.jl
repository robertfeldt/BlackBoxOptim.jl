function rosenbrock(x)
  return( sum( 100*( x[2:end] - x[1:end-1].^2 ).^2 + ( x[1:end-1] - 1 ).^2 ) )
end

function ellipsoid(x) 
  res = 0
  for(i in 1:length(x))
    res += sum(x[1:i])^2
  end
  res
end

function rastrigin(x)
  D = length(x)
  10 * D + sum( x.^2 ) - 10 * sum( cos( 2 * π * x ) )
end

function griewank(x)
  n = length(x)
  1 + (1/4000)*sum(x.^2) - prod(cos(x ./ sqrt(1:n)))
end

facts("compare_optimizers") do
  context("comparing optimizers on more than one problem") do
    ranks, fitnesses = BlackBoxOptim.compare_optimizers([
      (rosenbrock, (-5.0, 5.0)), 
      (rastrigin, (-5.12, 5.12)), 
      (ellipsoid, (-5.0, 5.0)),
      (griewank, (-600.0, 600.0))]; 
      dimensions = 30, max_time = 10.0)
    @fact size(ranks, 2) => 4
  end
end