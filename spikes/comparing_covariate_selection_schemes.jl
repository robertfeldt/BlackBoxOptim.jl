#
# Test a sparse rosenbrock, i.e. we have many factors but only a few are
# actually active in determining the response from the function.
#
function rosenbrock(x)
  n = length(x)
  return sum( 100*( x[2:n] - x[1:(n-1)].^2 ).^2 + ( x[1:(n-1)] - 1 ).^2 )
end
# FIXME this should be a more efficient implementation, but switching to it would reset all benchmarks
#rosenbrock(x) = sum(i -> 100*abs2(x[i+1] - x[i]^2) + abs2(x[i] - 1), Base.OneTo(length(x)-1))

function sphere(x)
  sum(x.^2)
end
# FIXME this should be a more efficient implementation, but switching to it would reset all benchmarks
#sphere(x) = sum(abs2, x)

mutable struct SparseShiftedProblem
  active_factors::Array{Int, 1}
  shift::Array{Float64, 2}
  func::Function

  function SparseShiftedProblem(n = 1000, k = 20, baseProblem = rosenbrock)
    afs = shuffle(collect(1:n))[1:k]
    new(afs, randn(k,1), baseProblem)
  end
end

function calc(sr::SparseShiftedProblem, x)
  sr.func(x[sr.active_factors] - sr.shift)
end

function sample_points_and_calc_response(n, s, func)
  xs = randn(n, s);
  y = zeros(s)
  [(y[i] = func(xs[:,i])) for i in 1:s]
  return xs, y
end

# Correlation selection is based on correlating each covariate with the response
# and keeping the covariates with the highest correlations with the response.
function correlation_selection(xs, y, maxSelected)
  correlations = y' * xs'
  find_top_largest_correlations = sort(correlations[:], rev=true)[1:maxSelected]
  index_of_largest_correlations = findall(in(find_top_largest_correlations), correlations)
  index_of_largest_correlations
end

# Absolute Correlation selection is based on correlating each covariate with the response
# and keeping the covariates with the highest absolute correlations with the response.
function abs_correlation_selection(xs, y, maxSelected)
  correlations = abs(y' * xs')
  find_top_largest_correlations = sort(correlations[:], rev=true)[1:maxSelected]
  index_of_largest_correlations = findall(in(find_top_largest_correlations), correlations)
  index_of_largest_correlations
end

# Stepwise filtering abs correlation selection. This is based on the idea that if
# there are complex interactive effects among the covariates the covariates that
# are part of an interaction should also have some individual effect. This is similar
# to the assumption of "effect heredity". Thus we do a stepwise filtering of the possible
# covariates to include by first eliminating half of them, then redoing the filtering
# on the interactions between the first selection, and so on.

# Random selection.
function random_selection(xs, y, maxSelected)
  n = size(xs, 1)
  shuffle(collect(1:n))[1:maxSelected]
end

# Call everything once to ensure they have been compiled.
a = rosenbrock([10.0, 4.5])
b = sphere([10.0, 4.5])
p = SparseShiftedProblem(100, 5, rosenbrock)
xs, y = sample_points_and_calc_response(100, 10, (x) -> calc(p, x))
s = random_selection(xs, y, 5)
s2 = correlation_selection(xs, y, 5)

# A problem is characterized mainly by:
#   1. How complex the interrelations are among the active covariates, basically determined by the objective func itself
#   2. Sparseness degree, i.e. how many of all available covariates that actually affects the response (0.0-1.0)
#   3. Dimensionality, the total number of covariates (100-10000)
#
# A selection scheme then saves an incomplete "trail" of evaluations of the objective function
# and uses them to predict which covariates are important. The size of this set is also
# important, as well as how it evolves over time.
#
# So we want to study how well a covariate subset selection (CSS) scheme predicts
# which covariates are important as we vary the four main parameters above.
# Dependent variables:
#   O1. Percentage of missed active covariates over NR repetitions
#   O2. Variation in missed covariates over NR repetitions
#   O3. Average time for a prediction
# Independent variables:
Problems = [rosenbrock, sphere]
Dimensionality = [1000, 10000, 100000]
ActualSparsenessDegree = [0.001, 0.01, 0.10, 0.25]
ArchiveSizeDegree = [0.1, 0.50, 1.00]
GuessedSparsenessDegree = [0.001, 0.01, 0.05, 0.10]

NumRepeats = 10
CovariateSelectionSchemes = [abs_correlation_selection, correlation_selection]
ncss = length(CovariateSelectionSchemes)

csvfile = open("comparing_covariate_selection_schemes.csv", "w")

using HypothesisTests

# Return the name of the best method of two methods given a Wilcoxon test.
# By default we are minimizing so better means smaller average, set minimizing
# to false if maximizing.
function best_of_if_better(x, y, names = ["1", "2"], minimizing = true, pvalue_treshold = 0.05)
  if pvalue(MannWhitneyUTest(x, y)) < pvalue_treshold
    op = minimizing ? < : >
    op(mean(x), mean(y)) ? names[1] : names[2]
  else
    "no_diff"
  end
end

println(csvfile, "CovariateSelectionScheme,Problem,Dimensionality,ActualSparseness,GuessedSparseness,ArchiveSize,MeanTime,StdMissed,MeanMissed")
for problem in Problems
  for dim in Dimensionality
    for asd in ActualSparsenessDegree
      for sizedegree in ArchiveSizeDegree
        for gsd in GuessedSparsenessDegree

          k = int(asd * dim)

          num_selected = int(gsd * dim)
          percent_missed = zeros(NumRepeats, ncss)
          times = zeros(NumRepeats, ncss)
          size = int(sizedegree * dim)

          for ir in 1:NumRepeats
            ssp = SparseShiftedProblem(dim, k, problem)
            xs, y = sample_points_and_calc_response(dim, size, (x) -> calc(ssp, x))

            for cssi in 1:ncss
              css = CovariateSelectionSchemes[cssi]
              tic()
              selected_factors = css(xs, y, num_selected)
              times[ir, cssi] = toq()
              missed = setdiff(Set(ssp.active_factors...), Set(selected_factors...))
              percent_missed[ir, cssi] = length(missed)/k*100.0
            end
          end

          for cssi in 1:ncss
            css = CovariateSelectionSchemes[cssi]
            println(csvfile, css, ",", problem, ",", dim, ",", asd, ",", gsd, ",", size, ",",
              mean(times[:,cssi]), ",",
              std(percent_missed[:,cssi]), ",", mean(percent_missed[:,cssi]))
            flush(stdout)
          end
          ms = mean(percent_missed, 1) # ms = mean(percent_missed, 1)
          best_method = best_of_if_better(percent_missed[:,1], percent_missed[:,2], CovariateSelectionSchemes)
          println(join([problem, dim, asd, gsd, size, minimum(ms), ms[2]-ms[1], best_method], ","))
        end
      end
    end
  end
end

# The conclusions from all these runs are:
#  For both problems:
#    1. Absolute correlation selection is almost always better than correlation selection (as expected)
#  For Sphere:
#    1. Only if the actual sparseness degree is very low (i.e. with < 30 actual active covariates) are the methods very good
#    2. And even so one must guess at least 10-50 times more active covariates
#  For Rosenbrock:
#    1. 1000: Only with an intermediate actual sparseness (0.01) are the methods very good
#    2. 10000: Only with very low actual sparseness (0.001) are the methods very good
#    3. 10000: Increasing the number of samples has quite a large effect i.e. lowers the missed percentage
#
# What is not clear is how good the method needs to be to actually give an advantage since they
# are applied in rounds. Maybe a small improvement can have quite a large effect in speeding up convergence.
#
# OTOH, we already have information about the correlation of different covariates in the diagonal
# of the covar matrix so maybe these additional selection algs are not needed...
