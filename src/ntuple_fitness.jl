"""
Base class for tuple-based fitness schemes.

Type parameters:
  * `N` is the number of the objectives
  * `F` is the type of each objective
  * `FA` is the actual type of the multi-objective fitness
  * `MIN` if objectives should be minimized or maximized
  * `AGG` the type of aggregator
"""
abstract type TupleFitnessScheme{N,F<:Number,FA,MIN,AGG} <: FitnessScheme{FA} end

@inline numobjectives(::TupleFitnessScheme{N}) where {N} = N
@inline fitness_eltype(::TupleFitnessScheme{N,F}) where {N,F} = F
@inline is_minimizing(::TupleFitnessScheme{N,F,FA,MIN}) where {N,F,FA,MIN} = MIN

@generated nafitness(::TupleFitnessScheme{N,F,NTuple{N,F}}) where {N,F} = ntuple(_ -> convert(F, NaN), N)
isnafitness(f::NTuple{N,F}, ::TupleFitnessScheme{N,F}) where {N,F} = any(isnan, f)

aggregate(f::NTuple{N,F}, fs::TupleFitnessScheme{N,F}) where {N,F} = fs.aggregator(f)

@inline is_better(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::TupleFitnessScheme{N,F,NTuple{N,F}}) where {N,F} = hat_compare(f1, f2, fs, -1) == -1
@inline is_worse(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::TupleFitnessScheme{N,F,NTuple{N,F}}) where {N,F} = hat_compare(f1, f2, fs, 1) == 1

"""
Pareto dominance for `N`-tuple (`N`≧1) fitnesses.

`aggregator::AGG` is a function mapping tuple fitness to a single numerical value.
Might be used for comparisons (or not, depending on the setup).
Always used when printing fitness vectors though.
"""
struct ParetoFitnessScheme{N,F<:Number,MIN,AGG} <: TupleFitnessScheme{N,F,NTuple{N,F},MIN,AGG}
    aggregator::AGG    # fitness aggregation function

    ParetoFitnessScheme{N,F,MIN,AGG}(; aggregator::AGG=sum) where {N, F<:Number, MIN, AGG} =
        new{N,F,MIN,AGG}(aggregator)

    ParetoFitnessScheme{N,F,MIN}(; aggregator::AGG=sum) where {N, F<:Number, MIN, AGG} =
        new{N,F,MIN,AGG}(aggregator)

    ParetoFitnessScheme{N,F}(; is_minimizing::Bool=true, aggregator::AGG=sum) where {N, F<:Number, AGG} =
        new{N,F,is_minimizing,AGG}(aggregator)

    ParetoFitnessScheme{N}(; fitness_type::Type{F} = Float64,
                            is_minimizing::Bool=true, aggregator::AGG=sum) where {N, F<:Number, AGG} =
        new{N,fitness_type,is_minimizing,AGG}(aggregator)
end

# comparison and for minimizing Pareto scheme
function hat_compare_pareto(u::NTuple{N,F}, v::NTuple{N,F}, expected::Int=0) where {N, F}
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

hat_compare(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::ParetoFitnessScheme{N,F,true}, expected::Int=0) where {N,F} =
    hat_compare_pareto(f1, f2, expected)
hat_compare(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::ParetoFitnessScheme{N,F,false}, expected::Int=0) where {N,F} =
    hat_compare_pareto(f2, f1, expected)

"""
ϵ-dominance for `N`-tuple (`N`≧1) fitnesses.

`aggregator::AGG` is a function mapping tuple fitness to a single numerical value.
Might be used for comparisons (or not, depending on the setup).
Always used when printing fitness vectors though.
"""
struct EpsDominanceFitnessScheme{N,F<:Number,MIN,AGG} <: FitnessScheme{NTuple{N,F}}
    ϵ::F              # ɛ-domination threshold
    aggregator::AGG    # fitness aggregation function

    function EpsDominanceFitnessScheme{N,F}(
            ϵ::F; is_minimizing::Bool=true, aggregator::AGG=sum) where {N, F<:Number, AGG}
        ϵ>0.0 || throw(ArgumentError("ϵ must be positive"))
        new{N,F,is_minimizing,AGG}(ϵ, aggregator)
    end
end

EpsDominanceFitnessScheme{N}(ϵ::F; fitness_type::Type{F} = Float64,
        is_minimizing::Bool=true, aggregator::AGG=sum) where {N, F<:Number, AGG} =
    EpsDominanceFitnessScheme{N,fitness_type}(ϵ; is_minimizing=is_minimizing, aggegator=aggregator)

# comparison for minimizing ϵ-dominance scheme
function hat_compare_ϵ(u::NTuple{N,F}, v::NTuple{N,F},
                       ϵ::F, expected::Int=0) where {N,F}
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

hat_compare(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::EpsDominanceFitnessScheme{N,F,true}, expected::Int=0) where {N,F} =
    hat_compare_ϵ(f1, f2, fs.ϵ, expected)
hat_compare(f1::NTuple{N,F}, f2::NTuple{N,F}, fs::EpsDominanceFitnessScheme{N,F,false}, expected::Int=0) where {N,F} =
    hat_compare_ϵ(f2, f1, fs.ϵ, expected)

# ϵ-index of the fitness component for minimizing scheme
@inline function ϵ_index(u::F, ϵ::F, ::Type{Val{true}}) where F
    if isnan(u)
        return (typemax(Int), zero(F))
    else
        u_div_ϵ = clamp(u/ϵ, convert(F, typemin(Int)), convert(F, typemax(Int)))
        ix = floor(Int, u_div_ϵ+10eps(F))
        return (ix, max(zero(F), u_div_ϵ-ix))
    end
end

# ϵ-index of the fitness component for maximizing scheme
@inline function ϵ_index(u::F, ϵ::F, ::Type{Val{false}}) where F
    if isnan(u)
        return (typemin(Int), zero(F))
    else
        u_div_ϵ = clamp(u/ϵ, convert(F, typemin(Int)), convert(F, typemax(Int)))
        ix = ceil(Int, u_div_ϵ)
        return (ix, ix-u_div_ϵ)
    end
end

# vectorized ϵ-index
@generated function ϵ_index(u::NTuple{N,F}, ϵ::Vector{F}, is_minimizing::Type{Val{MIN}}) where {N,F,MIN}
    quote
        pairs = Base.Cartesian.@ntuple $N i -> ϵ_index(u[i], ϵ[i], is_minimizing)
        ix = Base.Cartesian.@ntuple $N i -> pairs[i][1]
        sqrdist = zero(F)
        Base.Cartesian.@nexprs $N i -> sqrdist += pairs[i][2]^2
        return ix, sqrt(sqrdist)
    end
end

"""
ϵ-box indexed representation of the N-tuple fitness.

Used together with `EpsBoxDominanceFitnessScheme`.
"""
struct IndexedTupleFitness{N,F}
    orig::NTuple{N,F}       # original fitness
    agg::F                  # aggregated fitness
    index::NTuple{N,Int}    # ϵ-index vector
    dist::F                 # distance between ϵ-index vector and the original fitness

    function IndexedTupleFitness(u::NTuple{N,F}, agg::F, ϵ::Vector{F}, is_minimizing::Type{Val{MIN}}) where {N, F, MIN}
        ix, dist = ϵ_index(u, ϵ, is_minimizing)
        return new{N,F}(u, agg, ix, dist)
    end
end

IndexedTupleFitness(u::NTuple{N,F}, agg::F, ϵ::F, is_minimizing::Type{Val{MIN}}) where {N, F, MIN} =
    IndexedTupleFitness(u, agg, fill(ϵ, N), is_minimizing)

Base.convert(::Type{NTuple{N,F}}, fitness::IndexedTupleFitness{N,F}) where {N, F} = fitness.orig

@generated function nafitness(::Type{IndexedTupleFitness{N,F}}) where {N, F}
    quote
        IndexedTupleFitness(Base.Cartesian.@ntuple($N, _ -> convert($F, NaN)),
                            NaN, 1.0, Val{true})
    end
end

# comparison for minimizing ϵ-box dominance scheme
"""
Returns a tuple of `u` and `v` comparison:
  * `-1`: u≺v
  * `0`: u and v non-dominating
  * `1`: u≻v
and whether `u` index fully matches `v` index.
"""
function hat_compare_ϵ_box(
        u::IndexedTupleFitness{N,F},
        v::IndexedTupleFitness{N,F},
        is_minimizing::Bool=true, expected::Int=0) where {N,F}
    comp = 0
    @inbounds for (ui, vi) in zip(u.index, v.index)
        if ui > vi
            if comp == 0
                comp = 1
                if expected < 0 return (1, false) end
            elseif comp == -1
                return (0, false)  # non-dominated
            end
        elseif ui < vi
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
        return (comp, false)
    else
        uv_diff = u.dist - v.dist
        return (uv_diff < -10eps(F) ? -1 :
                uv_diff > 10eps(F) ? 1 : 0, true)
    end
end

function check_epsbox_ϵ(ϵ::Number, n::Int)
    ϵ>0.0 || throw(ArgumentError("ϵ must be positive"))
    return fill(ϵ, n)
end

function check_epsbox_ϵ(ϵ::AbstractVector{<:Number}, n::Int)
    length(ϵ)==n || throw(ArgumentError("The length of ϵ vector ($(length(ϵ))) does not match the specified fitness dimensions ($n)"))
    all(isposdef, ϵ) || throw(ArgumentError("ϵ must be positive"))
    return ϵ
end

"""
`EpsBoxDominanceFitnessScheme` defines ϵ-box dominance for
`N`-tuple (`N`≧1) fitnesses.
It operates with fitnesses of type `IndexedTupleFitness`.

`aggregator::AGG` is a function mapping tuple fitness to a single numerical value.
Might be used for comparisons (or not, depending on the setup).
Always used when printing fitness vectors though.
"""
struct EpsBoxDominanceFitnessScheme{N,F<:Number,MIN,AGG} <: TupleFitnessScheme{N,F,IndexedTupleFitness{N,F},MIN,AGG}
    ϵ::Vector{F}        # per-objective ɛ-domination thresholds
    aggregator::AGG     # fitness aggregation function

    EpsBoxDominanceFitnessScheme{N,F}(ϵ::Union{F,Vector{F}};
            is_minimizing::Bool=true, aggregator::AGG=sum) where {N,F<:Number,AGG} =
        new{N,F,is_minimizing,AGG}(check_epsbox_ϵ(ϵ, N), aggregator)

    EpsBoxDominanceFitnessScheme{N}(ϵ::Union{F,Vector{F}};
            is_minimizing::Bool=true, aggregator::AGG=sum) where {N,F<:Number,AGG} =
        new{N,F,is_minimizing,AGG}(check_epsbox_ϵ(ϵ, N), aggregator)
end

isnafitness(f::IndexedTupleFitness{N,F},
            fit_scheme::EpsBoxDominanceFitnessScheme{N,F}) where {N,F<:Number} =
    isnafitness(f.orig, fit_scheme)

EpsBoxDominanceFitnessScheme(fs::ParetoFitnessScheme{N,F}, ϵ::F=one(F)) where {N,F} =
    EpsBoxDominanceFitnessScheme{N,F}(ϵ, is_minimizing=is_minimizing(fs), aggregator=fs.aggregator)

EpsBoxDominanceFitnessScheme(fs::ParetoFitnessScheme{N,F}, ϵ::Vector{F}) where {N,F} =
    EpsBoxDominanceFitnessScheme{N,F}(ϵ, is_minimizing=is_minimizing(fs), aggregator=fs.aggregator)

EpsBoxDominanceFitnessScheme(fs::EpsDominanceFitnessScheme{N,F}, ϵ::Union{F,Vector{F}}=fs.ϵ) where {N,F} =
    EpsBoxDominanceFitnessScheme{N,F}(ϵ, is_minimizing=is_minimizing(fs), aggregator=fs.aggregator)

ParetoFitnessScheme(fs::EpsBoxDominanceFitnessScheme{N,F}) where {N,F} =
  ParetoFitnessScheme{N,F}(is_minimizing=is_minimizing(fs), aggregator=fs.aggregator)

Base.convert(::Type{IndexedTupleFitness{N,F}}, fitness::NTuple{N,F},
             fs::EpsBoxDominanceFitnessScheme{N,F}) where {N,F} =
    IndexedTupleFitness(fitness, aggregate(fitness, fs), fs.ϵ, Val{is_minimizing(fs)})

Base.convert(::Type{IndexedTupleFitness}, fitness::NTuple{N,F},
             fs::EpsBoxDominanceFitnessScheme{N,F}) where {N,F} =
    IndexedTupleFitness(fitness, aggregate(fitness, fs), fs.ϵ, Val{is_minimizing(fs)})

Base.convert(::Type{NTuple{N,F}}, fitness::IndexedTupleFitness{N,F},
             fs::EpsBoxDominanceFitnessScheme{N,F}) where {N,F} = fitness.orig

hat_compare(u::IndexedTupleFitness{N,F}, v::IndexedTupleFitness{N,F},
            fs::EpsBoxDominanceFitnessScheme{N,F}, expected::Int=0) where {N,F} =
    hat_compare_ϵ_box(u, v, is_minimizing(fs), expected)

hat_compare(u::NTuple{N,F}, v::IndexedTupleFitness{N,F},
            fs::EpsBoxDominanceFitnessScheme{N,F}, expected::Int=0) where {N,F} =
    hat_compare(convert(IndexedTupleFitness{N,F}, u, fs), v, fs, expected)

hat_compare(u::IndexedTupleFitness{N,F}, v::NTuple{N,F},
            fs::EpsBoxDominanceFitnessScheme{N,F}, expected::Int=0) where {N,F} =
    hat_compare(u, convert(IndexedTupleFitness{N,F}, v, fs), fs, expected)

hat_compare(u::NTuple{N,F}, v::NTuple{N,F},
            fs::EpsBoxDominanceFitnessScheme{N,F}, expected::Int=0) where {N,F} =
    hat_compare(convert(IndexedTupleFitness{N,F}, u, fs),
                convert(IndexedTupleFitness{N,F}, v, fs), fs, expected)

# special overload that strips index equality flag
(hc::HatCompare{FS})(u::IndexedTupleFitness{N,F},
                     v::IndexedTupleFitness{N,F}) where {FS<:EpsBoxDominanceFitnessScheme,N,F} =
    hat_compare(u, v, hc.fs)[1]
