using BlackBoxOptim

# From Hans W. Borchers in 
# https://groups.google.com/forum/#!topic/julia-opt/HltM-OGjueo
function trefethen(p)
  x = p[1]; y = p[2]
  return exp(sin(50 * x)) + sin(60 * exp(y)) + sin(70 * sin(x)) + 
             sin(sin(80 * y)) - sin(10 * (x + y)) + (x^2 + y^2)/4
end

p = BlackBoxOptim.fixeddim_problem(trefethen; range = (-1.0, 1.0), dims = 2, name = "trefethen")

# To just run one optimization:
bboptimize(p; max_time = 2.0, method = :generating_set_search, 
  parameters = {:SaveFitnessTraceToCsv => true})

# To compare several methods in several repetitions:
#BlackBoxOptim.repeated_bboptimize(5, p, 2, [
#  :generating_set_search, 
#  :adaptive_de_rand_1_bin_radiuslimited,
#  :random_search,
#  :separable_nes], 
#  5.0)