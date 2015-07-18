# A FitnessScheme is a specific way in which fitness vectors/values are
# aggregated, compared and presented. A fitness represents the score of one
# and the same individual on one or a set of evaluations. A scheme is a specific
# way in which these scores are considered in a coherent way.
# Type parameter F specifies the type of fitness values.
abstract FitnessScheme{F}

fitness_type{F}(::Type{FitnessScheme{F}}) = F
fitness_type{F}(::FitnessScheme{F}) = F
fitness_type{FS<:FitnessScheme}(::Type{FS}) = fitness_type(super(FS))
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
immutable ScalarFitnessScheme{MIN} <: RatioFitnessScheme{Float64}
end

const MinimizingFitnessScheme = ScalarFitnessScheme{true}()
const MaximizingFitnessScheme = ScalarFitnessScheme{false}()

is_minimizing{MIN}(::ScalarFitnessScheme{MIN}) = MIN
nafitness(::ScalarFitnessScheme) = NaN
isnafitness(f::Float64, ::ScalarFitnessScheme) = isnan(f)
numobjectives(::ScalarFitnessScheme) = 1

# Aggregation is just the identity function for scalar fitness
aggregate(fitness, ::ScalarFitnessScheme) = fitness

is_better(f1::Float64, f2::Float64, scheme::ScalarFitnessScheme{true}) = f1 < f2
is_better(f1::Float64, f2::Float64, scheme::ScalarFitnessScheme{false}) = f1 > f2

# Complex-valued fitness
# FIXME what is isbetter() for ComplexFitnessScheme
immutable ComplexFitnessScheme <: FitnessScheme{Complex128}
end

# FIXME do we need it? it might be confused with problem's fitness bounds 
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
immutable VectorFitnessScheme{MIN,N,AGG} <: RatioFitnessScheme{Vector{Float64}}
  # Function mapping a fitness array to a single numerical value. Might be used
  # for comparisons (or not, depending on setup). Always used when printing
  # fitness vectors though.
  aggregate::AGG
end

numobjectives{MIN,N}(::VectorFitnessScheme{MIN,N}) = N
is_minimizing{MIN,N}(fs::VectorFitnessScheme{MIN,N}) = MIN

nafitness{MIN,N}(::VectorFitnessScheme{MIN,N}) = fill(NaN, N)
isnafitness{MIN,N}(f::Vector{Float64}, ::VectorFitnessScheme{MIN,N}) = any(isnan(f)) # or all?

aggregate(fitness, fs::VectorFitnessScheme) = fs.aggregate(fitness)

# Fitness scheme that minimizes the sum of objectives
function vector_fitness_scheme_min(nobjectives::Int, aggregate = sum)
  VectorFitnessScheme{true, nobjectives, Function}(aggregate)
end

# Fitness scheme that maximizes the sum of objectives
function vector_fitness_scheme_max(nobjectives::Int, aggregate = sum)
  VectorFitnessScheme{false, nobjectives, Function}(aggregate)
end

# FIXME now it's here just to avoid undeclared types
type NewFitness end
