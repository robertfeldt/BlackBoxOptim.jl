# A FitnessScheme defines what constitutes a fitness, 
# how they are aggregated, and compared/ranked.
abstract FitnessScheme

# A FloatVectorFitness scheme has individual fitness scores (at least 1) in 
# an array.
type FloatVectorFitness <: FitnessScheme
  # Function mapping a fitness array to a single numerical value. Might be used
  # for comparisons (or not, depending on setup). Always used when printing
  # fitness vectors though.
  aggregate::Function

  # Worst aggregate fitness value.
  worst_fitness::FloatingPoint
end

hat_compare(a1::FloatingPoint, a2::FloatingPoint) = (a1 < a2) ? -1 : ((a1 > a2) ? 1 : 0)

# Hat comparison function that indicates which of fitness f1 and f2 is the better.
# Returns -1 if f1 is better than f2, 1 if f2 is better than f1 and
# 0 if they are equal.
function hat_compare(f1, f2, scheme::FitnessScheme)
  agg = scheme.aggregate
  a1, a2 = agg(f1), agg(f2)
  hat_compare(agg(f1), agg(f2))
end

isbetter(f1, f2, scheme::FitnessScheme) = hat_compare(f1, f2, scheme) == -1
isworse(f1, f2, scheme::FitnessScheme) = hat_compare(f1, f2, scheme) == 1
samefitness(f1, f2, scheme::FitnessScheme) = hat_compare(f1, f2, scheme) == 0

# For minimization we just pass the aggregator on.
function float_vector_scheme_min(aggregator = sum)
  FloatVectorFitness(aggregator, Inf)
end

# For maximization we need to set a different aggregator.
function float_vector_scheme_max(agg = ((fs) -> -1 * sum(fs)))
  FloatVectorFitness(agg, -Inf)
end
