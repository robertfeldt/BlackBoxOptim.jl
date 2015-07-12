# A FitnessScheme is a specific way in which fitness vectors/values are
# aggregated, compared and presented. A fitness represents the score of one
# and the same individual on one or a set of evaluations. A scheme is a specific
# way in which these scores are considered in a coherent way.
# Type parameter F specifies the type of fitness values.
abstract FitnessScheme{F}

fitness_type{F}(::Type{FitnessScheme{F}}) = F
fitness_type{F}(::FitnessScheme{F}) = F
#fitness_type{FS<:FitnessScheme}(::Type{FS}) = fitness_type(super(FS))
#fitness_type{FS<:FitnessScheme}(::FS) = fitness_type(FS)

if VERSION >= v"0.4.0-dev+1258" # FIXME remove version check once v0.4 is released
Base.call{F}(fs::FitnessScheme{F}, x::F, y::F) = is_better(x, y, fs)
fitness_scheme_lt(fs::FitnessScheme) = fs
else
fitness_scheme_lt(fs::FitnessScheme) = (x,y) -> is_better(x, y, fs)
end

# In a RatioFitnessScheme the fitness values can be ranked on a ratio scale so
# we need not rank them based on pairwise comparisons. The default scale used
# is the aggregate of the fitness values.
# FIXME
abstract RatioFitnessScheme{F} <: FitnessScheme{F}

# Fitness using a single floating value.
# The boolean type parameter specifies if smaller fitness is better (MIN=true)
# or worse (MIN=false).
immutable ScalarFitness{MIN} <: RatioFitnessScheme{Float64}
end

is_minimizing{MIN}(::ScalarFitness{MIN}) = MIN
nafitness(::ScalarFitness) = NaN
isnafitness(f::Float64, ::ScalarFitness) = isnan(f)
numobjectives(::ScalarFitness) = 1

# Aggregation is just the identity function for scalar fitness
aggregate(fitness, ::ScalarFitness) = fitness

is_better(f1::Float64, f2::Float64, scheme::ScalarFitness{true}) = f1 < f2
is_better(f1::Float64, f2::Float64, scheme::ScalarFitness{false}) = f1 > f2

# Complex-valued fitness
# FIXME what is isbetter() for ComplexFitness
immutable ComplexFitness <: FitnessScheme{Complex128}
end

worst_fitness(fs::FitnessScheme) = is_minimizing(fs) ? Inf : (-Inf)
best_fitness(fs::FitnessScheme) = -worst_fitness(fs)

hat_compare(a1::Number, a2::Number) = (a1 < a2) ? -1 : ((a1 > a2) ? 1 : 0)

# Hat comparison function that indicates which of fitness f1 and f2 is the better.
# Returns -1 if f1 is better than f2, 1 if f2 is better than f1 and
# 0 if they are equal.
function hat_compare(f1, f2, s::RatioFitnessScheme)
  if is_minimizing(s)
    hat_compare(aggregate(f1, s), aggregate(f2, s))
  else
    hat_compare(aggregate(f2, s), aggregate(f1, s))
  end
end

is_better(f1, f2, scheme::FitnessScheme) = hat_compare(f1, f2, scheme) == -1
is_worse(f1, f2, scheme::FitnessScheme) = hat_compare(f1, f2, scheme) == 1
same_fitness(f1, f2, scheme::FitnessScheme) = hat_compare(f1, f2, scheme) == 0

# All VectorFitness scheme has N individual fitness scores (at least 1) in
# an array and could be aggregated to a float value.
# FIXME
immutable VectorFitness{MIN,N,AGG} <: RatioFitnessScheme{Vector{Float64}}
  # Function mapping a fitness array to a single numerical value. Might be used
  # for comparisons (or not, depending on setup). Always used when printing
  # fitness vectors though.
  aggregate::AGG
end

numobjectives{MIN,N}(::VectorFitness{MIN,N}) = N
is_minimizing{MIN,N}(fs::VectorFitness{MIN,N}) = MIN

nafitness{MIN,N}(::VectorFitness{MIN,N}) = fill(NaN, N)
isnafitness{MIN,N}(f::Vector{Float64}, ::VectorFitness{MIN,N}) = any(isnan(f)) # or all?

aggregate(fitness, fs::VectorFitness) = fs.aggregate(fitness)

# Fitness scheme that minimizes the sum of objectives
function vector_fitness_scheme_min(nobjectives::Int, aggregate = sum)
  VectorFitness{true, nobjectives, Function}(aggregate)
end

# Fitness scheme that maximizes the sum of objectives
function vector_fitness_scheme_max(nobjectives::Int, aggregate = sum)
  VectorFitness{false, nobjectives, Function}(aggregate)
end

# FIXME now it's here just to avoid undeclared types
type NewFitness end
