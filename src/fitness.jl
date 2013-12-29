# A FitnessScheme is a specific way in which fitness vectors/values are
# aggregated, compared and presented. A fitness represents the score of one 
# and the same individual on one or a set of evaluations. A scheme is a specific 
# way in which these scores are considered in a coherent way.
abstract FitnessScheme

# In a RatioFitnessScheme the fitness values can be ranked on a ratio scale so
# we need not rank them based on pairwise comparisons. The default scale used
# is the aggregate of the fitness values.
abstract RatioFitnessScheme <: FitnessScheme

worst_fitness(fs::FitnessScheme) = is_minimizing(fs) ? Inf : (-Inf)
best_fitness(fs::FitnessScheme) = -worst_fitness(fs)
is_minimizing(fs::FitnessScheme) = true # Default is to minimize, override if not

hat_compare(a1::FloatingPoint, a2::FloatingPoint) = (a1 < a2) ? -1 : ((a1 > a2) ? 1 : 0)

# Hat comparison function that indicates which of fitness f1 and f2 is the better.
# Returns -1 if f1 is better than f2, 1 if f2 is better than f1 and
# 0 if they are equal.
function hat_compare(f1, f2, s::FitnessScheme)
  if is_minimizing(s)
    hat_compare(aggregate(f1, s), aggregate(f2, s))
  else
    hat_compare(aggregate(f2, s), aggregate(f1, s))
  end
end

is_better(f1, f2, scheme::FitnessScheme) = hat_compare(f1, f2, scheme) == -1
is_worse(f1, f2, scheme::FitnessScheme) = hat_compare(f1, f2, scheme) == 1
same_fitness(f1, f2, scheme::FitnessScheme) = hat_compare(f1, f2, scheme) == 0

# A FloatFitness scheme is simply a single floating value.
type FloatFitness <: RatioFitnessScheme
  minimize::Bool

  FloatFitness(minimize = true) = begin
    new(minimize)
  end
end

is_minimizing(fs::FloatFitness) = fs.minimize

# Aggregation is just the identity function for float fitness values.
aggregate(fitness, fs::FloatFitness) = fitness

# All FloatVectorFitness scheme has individual fitness scores (at least 1) in 
# an array.
type FloatVectorFitness <: RatioFitnessScheme
  # Function mapping a fitness array to a single numerical value. Might be used
  # for comparisons (or not, depending on setup). Always used when printing
  # fitness vectors though.
  aggregate::Function
end

aggregate(fitness, fs::FloatVectorFitness) = fs.aggregate(fitness)

# For minimization we just pass the aggregator on.
function float_vector_scheme_min(aggregator = sum)
  FloatVectorFitness(aggregator)
end

# For maximization we need to set a different aggregator.
function float_vector_scheme_max(agg = ((fs) -> -1 * sum(fs)))
  FloatVectorFitness(agg)
end