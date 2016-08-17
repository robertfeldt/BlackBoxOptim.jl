"""
  `FitnessScheme` defines how fitness vectors/values are
  compared, presented and aggregated.
  `Fitness` represents the score of one and the same individual
  on one or a set of evaluations. A scheme is a specific
  way to consider the scores in a coherent way.
  Type parameter `F` specifies the type of fitness values.

  `FitnessScheme` could also be used as a function that defines the fitness
  ordering, i.e. `fs(x, y) == true` iff fitness `x` is better than `y`.
"""
abstract FitnessScheme{F}

"""
  `fitness_type(fs::FitnessScheme)`
  `fitness_type(fs_type::Type{FitnessScheme})`

  Get the type of fitness values for fitness scheme `fs`.
"""
fitness_type{F}(::Type{FitnessScheme{F}}) = F
fitness_type{F}(::FitnessScheme{F}) = F
fitness_type{FS<:FitnessScheme}(::Type{FS}) = fitness_type(supertype(FS))
fitness_eltype{F<:Number}(::Type{FitnessScheme{F}}) = F
fitness_eltype{F<:Number}(::FitnessScheme{F}) = F
#fitness_type{FS<:FitnessScheme}(::FS) = fitness_type(FS)

# trivial convert() between calculated and archived fitness
Base.convert{F}(::Type{F}, fit::F, fit_scheme::FitnessScheme{F}) = fit

# ordering induced by the fitness scheme
# FIXME enable once v0.5 issue #14919 is fixed
# @compat (fs::FitnessScheme{F}){F}(x::F, y::F) = is_better(x, y, fs)

"""
  In `RatioFitnessScheme` the fitness values can be ranked on a ratio scale so
  pairwise comparisons are not required.
  The default scale used is the aggregate of the fitness components.
"""
# FIXME
abstract RatioFitnessScheme{F} <: FitnessScheme{F}

"""
  `Float64`-valued scalar fitness scheme.
  The boolean type parameter `MIN` specifies if smaller fitness values
  are better (`true`) or worse (`false`).
"""
immutable ScalarFitnessScheme{MIN} <: RatioFitnessScheme{Float64}
end

const MinimizingFitnessScheme = ScalarFitnessScheme{true}()
const MaximizingFitnessScheme = ScalarFitnessScheme{false}()

is_minimizing{MIN}(::ScalarFitnessScheme{MIN}) = MIN
nafitness{F<:Number}(::Type{F}) = convert(F, NaN)
@inline nafitness(fs::FitnessScheme) = nafitness(fitness_type(fs))
isnafitness{F<:Number}(f::F, ::RatioFitnessScheme{F}) = isnan(f)
numobjectives{F<:Number}(::RatioFitnessScheme{F}) = 1

""" Aggregation is just the identity function for scalar fitness. """
aggregate{F<:Number}(fitness::F, ::RatioFitnessScheme{F}) = fitness

is_better(f1::Float64, f2::Float64, scheme::ScalarFitnessScheme{true}) = f1 < f2
is_better(f1::Float64, f2::Float64, scheme::ScalarFitnessScheme{false}) = f1 > f2

""" Complex-valued fitness. """
# FIXME what is isbetter() for ComplexFitnessScheme
immutable ComplexFitnessScheme <: FitnessScheme{Complex128}
end

# FIXME do we need it? it might be confused with problem's fitness bounds
worst_fitness(fs::FitnessScheme) = is_minimizing(fs) ? Inf : (-Inf)
best_fitness(fs::FitnessScheme) = -worst_fitness(fs)

hat_compare(a1::Number, a2::Number) = (a1 < a2) ? -1 : ((a1 > a2) ? 1 : 0)

"""
  Check whether `f1` or `f2` fitness is better.

  Returns `-1` if `f1` is better than `f2`,
          `1` if `f2` is better than `f1` and
          `0` if they are equal.
"""
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

immutable HatCompare{FS<:FitnessScheme}
    fs::FS

    @compat (::Type{HatCompare}){FS<:FitnessScheme}(fs::FS) = new{FS}(fs)
end

@compat (hc::HatCompare){F}(x::F, y::F) = hat_compare(x, y, hc.fs)
