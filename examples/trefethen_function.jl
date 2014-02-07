using BlackBoxOptim

# From Hans W. Borchers in 
# https://groups.google.com/forum/#!topic/julia-opt/HltM-OGjueo
function ftrefethen(p)
  x = p[1]; y = p[2]
  return exp(sin(50 * x)) + sin(60 * exp(y)) + sin(70 * sin(x)) + 
             sin(sin(80 * y)) - sin(10 * (x + y)) + (x^2 + y^2)/4
end

bboptimize(ftrefethen, (-1.0, 1.0); dimensions = 2)