"""
An ordered set of dicts that are examined one after another to find the parameter value.
Returns nothing if no param setting is found in any of the dicts.
"""
mutable struct DictChain{K,V} <: AbstractDict{K,V}
    dicts::Vector{AbstractDict{K,V}}  # First dict is searched first and then in order until the last

    DictChain{K,V}(dicts::Vector{AbstractDict{K,V}}) where {K,V} = new{K,V}(dicts)

    # empty dicts vector
    DictChain{K,V}() where {K,V} = new(Vector{AbstractDict{K,V}}())

    DictChain{K,V}(dicts::AbstractDict{K,V}...) where {K,V} =
        new(AbstractDict{K,V}[dict for dict in dicts])
end

# (Associative{K,V}...) ctor is not triggered, so here's 1-, 2- and 3-argument
# versions of it, until there's a better way
DictChain(dict::AbstractDict{K,V}) where {K,V} = DictChain{K,V}(dict)
DictChain(dict1::AbstractDict{K,V}, dict2::AbstractDict{K,V}) where {K,V} =
    DictChain{K,V}(dict1, dict2)
DictChain(dict1::AbstractDict{K,V},
          dict2::AbstractDict{K,V},
          dict3::AbstractDict{K,V}) where {K,V} =
    DictChain{K,V}(dict1, dict2, dict3)

function Base.show(io::IO, dc::DictChain)
    print(io, typeof(dc), "[")
    for (i, dict) in enumerate(dc.dicts)
        (i > 1) && print(io, ",")
        show(io, dict)
    end
    print(io, "]")
end

Base.show(io::IO, ::MIME{Symbol("text/plain")}, dc::DictChain) = show(io, dc)

function Base.getindex(p::DictChain, key)
    for d in p.dicts
        haskey(d, key) && return getindex(d, key)
    end
    throw(KeyError(key))
end

function Base.setindex!(p::DictChain{K,V}, value, key) where {K,V}
    isempty(p.dicts) && push!(p.dicts, Dict{K,V}()) # add new dictionary on top
    setindex!(first(p.dicts), value, key)
end

function Base.haskey(p::DictChain, key)
    for d in p.dicts
        haskey(d, key) && return true
    end
    return false
end

Base.get(p::DictChain, key, default = nothing) =
    return haskey(p, key) ? getindex(p, key) : default

# In a merge the last parameter should be prioritized since this is the way
# the normal Julia merge() of dicts works.
Base.merge(p1::DictChain{K,V}, p2::Dict{K,V}) where {K,V} = DictChain{K,V}([p2; p1.dicts])
Base.merge(p1::DictChain{K,V}, p2::DictChain{K,V}) where {K,V} = DictChain{K,V}([p2.dicts; p1.dicts])
Base.merge(p1::Dict{K,V}, p2::DictChain{K,V}) where {K,V} = DictChain{K,V}([p2.dicts; p1])

function Base.merge!(dc::DictChain{K,V}, d::Dict{K,V}) where {K,V}
    insert!(dc.dicts, 1, d)
    return dc
end

function Base.merge!(d::Dict{K,V}, dc::DictChain{K,V}) where {K,V}
    for dict in reverse(dc.dicts)
        merge!(d, dict)
    end
    return d
end

# since Base.merge(Dict, Dict) generates Dict, we need another name
# for operation that generates DictChain.
# The difference between chain() and merge() for other argument types is that
# merge() grows horizontally (extends dicts vector),
# whereas chain() grows vertically
# (references its two arguments in the new 2-element dicts vector)
chain(p1::AbstractDict{K,V}, p2::AbstractDict{K,V}) where {K,V} = DictChain(p2, p1)
chain(p1::AbstractDict{K,V}, p2::AbstractDict{K,V}...) where {K,V} = DictChain(chain(p2...), p1)

flatten(d::AbstractDict) = d
flatten(d::DictChain{K,V}) where {K,V} = convert(Dict{K,V}, d)

Base.convert(::Type{Dict{K,V}}, dc::DictChain{K,V}) where {K,V} = merge!(Dict{K,V}(), dc)

function Base.delete!(p::DictChain, key)
    for d in p.dicts
        delete!(d, key)
    end
    return p
end

"""
The parameters storage type for `BlackBoxOptim`.
"""
const Parameters = AbstractDict{Symbol,Any}

"""
The default parameters storage in `BlackBoxOptim`.
"""
const ParamsDict = Dict{Symbol,Any}
const ParamsDictChain = DictChain{Symbol,Any}

"""
The default placeholder value for parameters argument.
"""
const EMPTY_PARAMS = ParamsDict()

kwargs2dict(kwargs::Iterators.Pairs{Symbol}) = ParamsDict(kwargs)
