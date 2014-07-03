using BlackBoxOptim

# Another problem from Hans W. Borchers, constrained 15-points:
function points15(x)
    d = 2.0
    for i = 1:14, j = (i+1):15
        s = (x[i]-x[j])^2 + (x[i+15]-x[j+15])^2 + (x[i+30]-x[j+30])^2
        s = sqrt(s)
        if s < d
          d = s
        end
    end
    return -d
end

Dims = 15 * 3

# The best known minimum is around -0.62 but we do not clearly know if there are better
# ones so we aim somewhat lower...
p = BlackBoxOptim.fixeddim_problem(points15; range = (0.0, 1.0), dims = Dims, name = "constrained points 15", fmins = [-0.64])

# This is a hard problem so run for some time...
MaxMinutes = 0.50

# Run the AdaptiveDE once:
#bboptimize(points15, (0.0, 1.0); dimensions = 45, max_time = MaxMinutes * 60)

# Run several methods several times:
BlackBoxOptim.repeated_bboptimize(2, p, Dims, [
  :generating_set_search, 
  #:adaptive_de_rand_1_bin_radiuslimited,
  #:random_search,
  :separable_nes], 
  MaxMinutes * 60)

#compare_optimizers(points15; search_range = (0.0, 1.0), dimensions = 45, 
#  max_time = MaxMinutes * 60, 
#  methods = [
#    :generating_set_search,
#    :adaptive_de_rand_1_bin_radiuslimited,
#  ]
#)

# To run XNES on it we need to add a penalty for going outside the (0,1) box.
function penalty(x, range)
  lower = x .< range[1]
  higher = x .> range[2]
  100 * (norm(range[1] - x[lower]) + norm(x[higher] - range[2]))
end

penalized_points15(x) = points15(x) + penalty(x, (0.0, 1.0))

# Run XNES
#best, fitness = bboptimize(penalized_points15, (0.0, 1.0); dimensions = 45, 
#  max_time = MaxMinutes * 60, method = :xnes)

#println("points15 fitness = ", points15(best))
