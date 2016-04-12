"""
  Base class for tuple-based fitness schemes.
"""
abstract TupleFitnessScheme{N,F<:Number,MIN,AGG} <: FitnessScheme{NTuple{N,F}}

@inline numobjectives{N}(::TupleFitnessScheme{N}) = N
@inline fitness_eltype{N,F}(::TupleFitnessScheme{N,F}) = F
@inline is_minimizing{N,F,MIN}(::TupleFitnessScheme{N,F,MIN}) = MIN

@generated nafitness{N,F}(::TupleFitnessScheme{N,F}) = ntuple(_ -> convert(F, NaN), Val{N})
isnafitness{N,F}(f::NTuple{N,F}, ::TupleFitnessScheme{N,F}) = any(isnan, f) # or any?

aggregate{N,F}(f::NTuple{N,F}, fs::TupleFitnessScheme{N,F}) = fs.aggregator(f)

@inline is_better{N,F}(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::TupleFitnessScheme{N,F}) = hat_compare(f1, f2, fs, -1) == -1
@inline is_worse{N,F}(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::TupleFitnessScheme{N,F}) = hat_compare(f1, f2, fs, 1) == 1

"""
  Pareto dominance for `N`-tuple (`N`≧1) fitnesses.

  `aggregator::AGG` is a function mapping tuple fitness to a single numerical value.
  Might be used for comparisons (or not, depending on the setup).
  Always used when printing fitness vectors though.
"""
immutable ParetoFitnessScheme{N,F<:Number,MIN,AGG} <: TupleFitnessScheme{N,F,MIN,AGG}
    aggregator::AGG    # fitness aggregation function

    Base.call{N,F<:Number,AGG}(::Type{ParetoFitnessScheme{N,F}};
                                is_minimizing::Bool=true, aggregator::AGG=sum) =
        new{N,F,is_minimizing,AGG}(aggregator)

    Base.call{N,F<:Number,AGG}(::Type{ParetoFitnessScheme{N}}; fitness_type::Type{F} = Float64,
                                is_minimizing::Bool=true, aggregator::AGG=sum) =
        new{N,fitness_type,is_minimizing,AGG}(aggregator)
end

# comparison and for minimizing Pareto scheme
function hat_compare_pareto{N,F}(u::NTuple{N,F}, v::NTuple{N,F}, expected::Int=0)
    res = 0
    @inbounds for i in 1:N
        delta = u[i] - v[i]
        if delta > 0.0
            if res == 0
                res = 1
                if expected == -1 return res end
            elseif res == -1
                return 0 # non-dominated
            end
        elseif delta < 0.0
            if res == 0
                res = -1
                if expected == 1 return res end
            elseif res == 1
                return 0 # non-dominated
            end
        end
    end
    return res
end

hat_compare{N,F}(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::ParetoFitnessScheme{N,F,true}, expected::Int=0) =
    hat_compare_pareto(f1, f2, expected)
hat_compare{N,F}(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::ParetoFitnessScheme{N,F,false}, expected::Int=0) =
    hat_compare_pareto(f2, f1, expected)

"""
  ϵ-dominance for `N`-tuple (`N`≧1) fitnesses.

  `aggregator::AGG` is a function mapping tuple fitness to a single numerical value.
  Might be used for comparisons (or not, depending on the setup).
  Always used when printing fitness vectors though.
"""
immutable EpsDominanceFitnessScheme{N,F<:Number,MIN,AGG} <: FitnessScheme{NTuple{N,F}}
    ϵ::F              # ɛ-domination threshold
    aggregator::AGG    # fitness aggregation function

    function Base.call{N,F<:Number,AGG}(::Type{EpsDominanceFitnessScheme{N,F}},
                                ϵ::F; is_minimizing::Bool=true, aggregator::AGG=sum)
        ϵ>0.0 || throw(ArgumentError("ϵ must be positive"))
        new{N,F,is_minimizing,AGG}(ϵ, aggregator)
    end

    Base.call{N,F<:Number,AGG}(::Type{EpsDominanceFitnessScheme{N}}, ϵ::F; fitness_type::Type{F} = Float64,
                               is_minimizing::Bool=true, aggregator::AGG=sum) =
        EpsDominanceFitnessScheme{N,fitness_type}(ϵ; is_minimizing=is_minimizing, aggegator=aggregator)
end

# comparison for minimizing ϵ-dominance scheme
function hat_compare_ϵ{N,F}(u::NTuple{N,F}, v::NTuple{N,F}, ϵ::F, expected::Int=0)
    res = 0 # true if any u[i] < v[i] + ϵ
    @inbounds for i in 1:N
        delta = u[i] - v[i] - ϵ
        if delta > 0.0
            if res == 0
                res = 1
                if expected == -1 return 1 end
            elseif res == -1
                return 0 # non-dominated
            end
        elseif delta < 0.0
            if res == 0
                res = -1
                if expected == 1 return -1 end
            elseif res == 1
                return 0 # non-dominated
            end
        end
    end
    return res
end

hat_compare{N,F}(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::EpsDominanceFitnessScheme{N,F,true}, expected::Int=0) =
    hat_compare_ϵ(f1, f2, fs.ϵ, expected)
hat_compare{N,F}(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::EpsDominanceFitnessScheme{N,F,false}, expected::Int=0) =
    hat_compare_ϵ(f2, f1, fs.ϵ, expected)

# ϵ-index of the fitness component for minimizing scheme
@inline function ϵ_index{F}(u::F, ϵ::F, ::Type{Val{true}})
    if isnan(u)
        return (typemax(Int), zero(F))
    else
        u_div_ϵ = clamp(u/ϵ, convert(F, typemin(Int)), convert(F, typemax(Int)))
        ix = floor(Int, u_div_ϵ+10eps(F))
        return (ix, max(zero(F), u_div_ϵ-ix))
    end
end

# ϵ-index of the fitness component for maximizing scheme
@inline function ϵ_index{F}(u::F, ϵ::F, ::Type{Val{false}})
    if isnan(u)
        return (typemin(Int), zero(F))
    else
        u_div_ϵ = clamp(u/ϵ, convert(F, typemin(Int)), convert(F, typemax(Int)))
        ix = ceil(Int, u_div_ϵ)
        return (ix, ix-u_div_ϵ)
    end
end

"""
    ϵ-box indexed representation of the N-tuple fitness.

    Used together with `EpsBoxDominanceFitnessScheme`.
"""
immutable IndexedTupleFitness{N,F}
    orig::NTuple{N,F}       # original fitness
    agg::F                  # aggregated fitness
    index::NTuple{N,Int}    # ϵ-index vector
    dist::F                 # distance between ϵ-index vector and the original fitness

    # TODO use @generated to avoid loops for N<=8?
    function Base.call{N,F,MIN}(::Type{IndexedTupleFitness}, u::NTuple{N,F}, agg::F, ϵ::Vector{F}, is_minimizing::Type{Val{MIN}})
        pairs = ntuple(i -> ϵ_index(u[i], ϵ[i], is_minimizing), Val{N}) # FIXME faster?
        return new{N,F}(u, agg, ntuple(i -> pairs[i][1], Val{N}),
                        sqrt(sum(pair::Tuple{Int,F} -> pair[2]^2, pairs)))
    end
    Base.call{N,F,MIN}(::Type{IndexedTupleFitness}, u::NTuple{N,F}, agg::F, ϵ::F, is_minimizing::Type{Val{MIN}}) =
        IndexedTupleFitness(u, agg, fill(ϵ, N), is_minimizing)
end

Base.convert{N,F}(::Type{NTuple{N,F}}, fitness::IndexedTupleFitness{N,F}) = fitness.orig

@generated nafitness{N,F}(::Type{IndexedTupleFitness{N,F}}) = IndexedTupleFitness(ntuple(_ -> convert(F, NaN), Val{N}), NaN, 1.0, Val{true})

# comparison for minimizing ϵ-box dominance scheme
"""
    Returns a tuple of `u` and `v` comparison:
      * `-1`: u≺v
      * `0`: u and v non-dominating
      * `1`: u≻v
    and whether `u` index fully matches `v` index.
"""
function hat_compare_ϵ_box{N,F}(u::IndexedTupleFitness{N,F}, v::IndexedTupleFitness{N,F}, is_minimizing::Bool=true, expected::Int=0)
    comp = 0
    @inbounds for i in 1:N
        if u.index[i] > v.index[i]
            if comp == 0
                comp = 1
                if expected < 0 return (1, false) end
            elseif comp == -1
                return (0, false)  # non-dominated
            end
        elseif u.index[i] < v.index[i]
            if comp == 0
                comp = -1
                if expected > 0 return (-1, false) end
            elseif comp == 1
                return (0, false) # non-dominated
            end
        end
    end
    if !is_minimizing
        comp = -comp
    end
    if comp != 0
        return (comp, comp==0)
    else
        uv_diff = u.dist - v.dist
        return (uv_diff < -10eps(F) ? -1 : uv_diff > 10eps(F) ? 1 : 0, comp==0)
    end
end

function check_epsbox_ϵ(ϵ::Number, n::Int)
    ϵ>0.0 || throw(ArgumentError("ϵ must be positive"))
    fill(ϵ, n)
end

function check_epsbox_ϵ{F<:Number}(ϵ::Vector{F}, n::Int)
    length(ϵ)==n || throw(ArgumentError("The length of ϵ vector ($(length(ϵ))) does not match the specified fitness dimensions ($n)"))
    all(isposdef, ϵ) || throw(ArgumentError("ϵ must be positive"))
    ϵ
end

"""
  `EpsBoxDominanceFitnessScheme` defines ϵ-box dominance for
  `N`-tuple (`N`≧1) fitnesses.
  It operates with fitnesses of type `IndexedTupleFitness`.

  `aggregator::AGG` is a function mapping tuple fitness to a single numerical value.
  Might be used for comparisons (or not, depending on the setup).
  Always used when printing fitness vectors though.
"""
immutable EpsBoxDominanceFitnessScheme{N,F<:Number,MIN,AGG} <: TupleFitnessScheme{N,F,MIN,AGG}
    ϵ::Vector{F}        # per-objective ɛ-domination thresholds
    aggregator::AGG     # fitness aggregation function

    Base.call{N,F<:Number,AGG}(::Type{EpsBoxDominanceFitnessScheme{N,F}}, ϵ::Union{F,Vector{F}};
                               is_minimizing::Bool=true, aggregator::AGG=sum) =
        new{N,F,is_minimizing,AGG}(check_epsbox_ϵ(ϵ, N), aggregator)

    Base.call{N,F<:Number,AGG}(::Type{EpsBoxDominanceFitnessScheme{N}}, ϵ::Union{F,Vector{F}};
                               is_minimizing::Bool=true, aggregator::AGG=sum) =
        new{N,F,is_minimizing,AGG}(check_epsbox_ϵ(ϵ, N), aggregator)
end

# overloads of the default behaviour, because the actual fitness type is not NTuple{N,F}
fitness_type{N,F}(::Type{EpsBoxDominanceFitnessScheme{N,F}}) = IndexedTupleFitness{N,F}
nafitness{N,F}(::EpsDominanceFitnessScheme{N,F}) = nafitness(IndexedTupleFitness{N,F})
isnafitness{N,F}(f::IndexedTupleFitness{N,F}, ::EpsBoxDominanceFitnessScheme{N,F}) = any(isnan, f.orig) # or any?

Base.convert{N,F,MIN}(::Type{EpsBoxDominanceFitnessScheme}, fs::ParetoFitnessScheme{N,F,MIN}, ϵ::F=one(F)) =
  EpsBoxDominanceFitnessScheme{N,F}(ϵ, is_minimizing=MIN, aggregator=fs.aggregator)

Base.convert{N,F,MIN}(::Type{EpsBoxDominanceFitnessScheme}, fs::ParetoFitnessScheme{N,F,MIN}, ϵ::Vector{F}) =
  EpsBoxDominanceFitnessScheme{N,F}(ϵ, is_minimizing=MIN, aggregator=fs.aggregator)

Base.convert{N,F,MIN}(::Type{EpsBoxDominanceFitnessScheme}, fs::EpsDominanceFitnessScheme{N,F,MIN}, ϵ::Union{F,Vector{F}}=fs.ϵ) =
  EpsBoxDominanceFitnessScheme{N,F}(ϵ, is_minimizing=MIN, aggregator=fs.aggregator)

Base.convert{N,F,MIN}(::Type{IndexedTupleFitness{N,F}}, fitness::NTuple{N,F},
                      fs::EpsBoxDominanceFitnessScheme{N,F,MIN}) =
    IndexedTupleFitness(fitness, aggregate(fitness, fs), fs.ϵ, Val{MIN})

Base.convert{N,F,MIN}(::Type{IndexedTupleFitness}, fitness::NTuple{N,F},
                      fs::EpsBoxDominanceFitnessScheme{N,F,MIN}) =
    IndexedTupleFitness(fitness, aggregate(fitness, fs), fs.ϵ, Val{MIN})

hat_compare{N,F,MIN}(u::IndexedTupleFitness{N,F}, v::IndexedTupleFitness{N,F},
                 fs::EpsBoxDominanceFitnessScheme{N,F,MIN}, expected::Int=0) =
    hat_compare_ϵ_box(u, v, MIN, expected)

hat_compare{N,F}(u::NTuple{N,F}, v::IndexedTupleFitness{N,F},
                 fs::EpsBoxDominanceFitnessScheme{N,F}, expected::Int=0) =
    hat_compare(convert(IndexedTupleFitness, u, fs), v, fs, expected)

hat_compare{N,F}(u::IndexedTupleFitness{N,F}, v::NTuple{N,F},
                 fs::EpsBoxDominanceFitnessScheme{N,F}, expected::Int=0) =
    hat_compare(u, convert(IndexedTupleFitness, v, fs), fs, expected)

hat_compare{N,F}(u::NTuple{N,F}, v::NTuple{N,F},
                 fs::EpsBoxDominanceFitnessScheme{N,F}, expected::Int=0) =
    hat_compare(convert(IndexedTupleFitness, u, fs),
                convert(IndexedTupleFitness, v, fs), fs, expected)

# special overload that strips index equality flag
Base.call{FS<:EpsBoxDominanceFitnessScheme,N,F}(hc::HatCompare{FS},
        u::IndexedTupleFitness{N,F}, v::IndexedTupleFitness{N,F}) = hat_compare(u, v, hc.fs)[1]
