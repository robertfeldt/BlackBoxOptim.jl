# Find optimal bin width for a histogram using the method described in:
#  Shimazaki H. and Shinomoto S., A method for selecting the bin size of a
#  time histogram Neural Computation (2007) Vol. 19(6), 1503-1527.

# We need the biased variance and not the unbiased which is implemented in Julia.
function biased_var(x)
  n = length(x)
  var(x) * (n - 1) / n
end

# Assumes the data is sorted.
function count_events_per_bin(delta, sorteddata)
  minv = sorteddata[1]
  numbins = floor(Int, (sorteddata[end] - minv) / delta) + 1
  counts = zeros(Int, numbins)

  current_bin = 1
  next_bin_start = minv + delta

  for data in sorteddata
    if data >= next_bin_start
      next_bin_start += delta
      current_bin += 1
    end
    counts[current_bin] += 1
  end

  counts
end

# The function to be optimized is a function taking the bin-width, delta, as
# the input, but it also needs the data to calc C(delta).
function c_objective(delta, data)
  bincounts = count_events_per_bin(delta, data)
  k = mean(bincounts)
  v = biased_var(bincounts)
  (2*k - v) / delta
end

using BlackBoxOptim

function find_optimal_bin_width(data; minwidth = nothing, maxwidth = nothing)
  sorted = sort(data)
  minvalue = sorted[1]
  maxvalue = sorted[end]
  range = maxvalue - minvalue

  min_width = (minwidth == nothing) ? range/200 : minwidth
  max_width = (maxwidth == nothing) ? range/2   : maxwidth

  best, fitness = bboptimize((delta) -> c_objective(delta[1], sorted), (min_width, max_width);
    dimensions = 1, iterations = 1e4)

  best[1], count_events_per_bin(best[1], sorted)
end
