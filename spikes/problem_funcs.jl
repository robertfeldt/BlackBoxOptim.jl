function rosenbrock(x)
  n = length(x)
  sum( 100*( x[2:n] - x[1:(n-1)].^2 ).^2 + ( x[1:(n-1)] - 1 ).^2 )
end

function sphere(x)
  sum(x .^ 2)
end

function cigar(x)
  x[1]^2 + 1e6 * sum(x[2:end].^2)
end

function cigtab(x)
  x[1]^2 + 1e8 * x[end]^2 + 1e4 * sum(x[2:(end-1)].^2)
end

function xtransform(n, f)
  move = randn(n, 1)
  transformed_f(x) = f(x .- move)
end