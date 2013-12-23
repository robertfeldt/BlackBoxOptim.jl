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

# From FEP paper:

# f1 is sphere
f1 = sphere

# f2 is schwefel2_22
function f2(x)
  absx = abs(x)
  sum(absx) + prod(absx)
end

# f3 in FEP paper is schwefel1_2
function schwefel1_2(x)
  s = 0
  for i in 1:length(x)
    partsum = sum(x[1:i])
    s+= (partsum^2)
  end
  s
end
f3 = schwefel1_2

# f4 is schwefel2_21
function schwefel2_21(x)
  maximum(abs(x))
end
f4 = schwefel2_21

# f5 is generalized rosenbrock

# f6 step function
function f6(x)
  sum(ceil(x + 0.5))
end

# f7 is quartic function with noise

# f8 is Generalized Schwefel 2_26

# f9 is GEneralized Rastrigin

# f10 Ackley's function

# f11 is Generalized Griewank

# 
function f8(x)
  sum( (-x) .* sin(sqrt(abs(x))) )
end

# This is a generator for the family of deceptive functions from the 
# Cuccu2011 paper on novelty-based restarts. We have vectorized it to allow
# more than 1D versions.
function deceptive_cuccu2011(l, w)
  f(x) = begin
    absx = abs(x)
    sumabsx = sum(absx)
    if sumabsx <= 1
      return sum(x^2)
    elseif sumabsx >= l+1
      return sum((absx - l)^2)
    else
      return sum(1 - 0.5 * (sin( (π * w * (absx - 1)) / l ))^2)
    end
  end
  return f
end

# As described in the Cuccu2011 paper.
function rastrigin(x)
  d = length(x)
  10*d + sum( x.^2 - 10 * cos( 2 * π * x ) )
end



function xtransform(n, f)
  move = randn(n, 1)
  transformed_f(x) = f(x .- move)
end