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
prob = BlackBoxOptim.minimization_problem(points15, "Constrained points 15",
        (0.0, 1.0), 45, -0.64)

# This is a hard problem so run for some time...
MaxMinutes = 0.50

# To run XNES et al on it we need to add a penalty for going outside the (0,1) box.
function penalty(x, range)
  lower = x .< range[1]
  higher = x .> range[2]
  100 * (norm(range[1] - x[lower]) + norm(x[higher] - range[2]))
end

penalized_points15(x) = points15(x) + penalty(x, (0.0, 1.0))

# Run XNES
result = bboptimize(penalized_points15; SearchRange = (0.0, 1.0), NumDimensions = Dims,
  MaxTime = MaxMinutes * 60, Method = :xnes)
best_xnes = best_candidate(result)

# Run DXNES
result = bboptimize(penalized_points15; SearchRange = (0.0, 1.0), NumDimensions = Dims,
  MaxTime = MaxMinutes * 60, Method = :dxnes)
best_dxnes = best_candidate(result)

# Run the default (Adaptive DE) method for same amount of time on penalized function
result = bboptimize(penalized_points15; SearchRange = (0.0, 1.0), NumDimensions = Dims,
  MaxTime = MaxMinutes * 60)
best_de_pen = best_candidate(result)

# But DE need no penalization since it respects search range bounds:
result = bboptimize(points15; SearchRange = (0.0, 1.0), NumDimensions = Dims,
  MaxTime = MaxMinutes * 60)
best_de = best_candidate(result)

# Run several methods several times:
BlackBoxOptim.repeated_bboptimize(2, prob, Dims, [
  :generating_set_search,
  :adaptive_de_rand_1_bin_radiuslimited,
  :random_search,
  ],
  MaxMinutes * 60)

println("xNES points15 fitness = ", points15(best_xnes))
println("dxNES points15 fitness = ", points15(best_dxnes))
println("Adaptive DE (default, penalized) points15 fitness = ", points15(best_de_pen))
println("Adaptive DE (default) points15 fitness = ", points15(best_de))
