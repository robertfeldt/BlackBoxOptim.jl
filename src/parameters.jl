# An ordered set of dicts which we traverse to find the parameter value.
# Returns nothing if no param setting is found in any of the dicts.
type DictChain{K,V} <: Associative{K,V}
  dicts::Vector{Associative{K,V}}  # First dict is searched first and then in order until the last

  DictChain(dicts::Vector{Associative{K,V}}) = new(dicts)

  # empty dicts vector
  DictChain() = new(Array(Associative{K,V}, 0))

  DictChain(dicts::Associative{K,V}...) = new(Associative{K,V}[dict for dict in dicts])
end

# (Associative{K,V}...) ctor is not triggered, so here's 1-, 2- and 3-argument
# versions of it, until there's a better way
DictChain{K,V}(dict::Associative{K,V}) = DictChain{K,V}(dict)
DictChain{K,V}(dict1::Associative{K,V}, dict2::Associative{K,V}) = DictChain{K,V}(dict1, dict2)
DictChain{K,V}(dict1::Associative{K,V}, dict2::Associative{K,V}, dict3::Associative{K,V}) = DictChain{K,V}(dict1, dict2, dict3)

function Base.show{K,V}(io::IO, dc::DictChain{K,V})
  print(io, typeof(dc), "[")
  for (i, dict) in enumerate(dc.dicts)
    if i > 1
      print(",")
    end
    show(io, dict)
  end
  print(io, "]")
end

function Base.getindex{K,V}(p::DictChain{K,V}, key::K)
  for d in p.dicts
    if haskey(d, key)
      return getindex(d, key)
    end
  end
  throw(KeyError(key))
end

# FIXME no type declaration for value due to 0.3 compatibility
function Base.setindex!{K,V}(p::DictChain{K,V}, value, key::K)
  if isempty(p.dicts)
    push!(p.dicts, Dict{K,V}()) # add new dictionary on top
  end
  setindex!(first(p.dicts), value, key)
end

import Base.haskey, Base.get, Base.delete!

function Base.haskey{K,V}(p::DictChain{K,V}, key::K)
  for d in p.dicts
    if haskey(d, key)
      return true
    end
  end
  return false
end

# FIXME no type declaration for default due to 0.3 compatibility
function get{K,V}(p::DictChain{K,V}, key::K, default = nothing)
  return haskey(p, key) ? getindex(p, key) : default
end

# In a merge the last parameter should be prioritized since this is the way
# the normal Julia merge() of dicts works.
Base.merge{K,V}(p1::DictChain{K,V}, p2::Dict{K,V}) = DictChain{K,V}([p2; p1.dicts])
Base.merge{K,V}(p1::DictChain{K,V}, p2::DictChain{K,V}) = DictChain{K,V}([p2.dicts; p1.dicts])
Base.merge{K,V}(p1::Dict{K,V}, p2::DictChain{K,V}) = DictChain{K,V}([p2.dicts; p1])

function Base.merge!{K,V}(dc::DictChain{K,V}, d::Dict{K,V})
  insert!(dc.dicts, 1, d); dc
end

function Base.merge!{K,V}(d::Dict{K,V}, dc::DictChain{K,V})
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
chain{K,V}(p1::Associative{K,V}, p2::Associative{K,V}) = DictChain(p2, p1)
chain{K,V}(p1::Associative{K,V}, p2::Associative{K,V}...) = DictChain(chain(p2...), p1)

flatten(d::Associative) = d
flatten{K,V}(d::DictChain{K,V}) = convert(Dict{K,V}, d)

# converts DictChain into Dict
Base.convert{K,V}(::Type{Dict{K,V}}, dc::DictChain{K,V}) = merge!(Dict{K,V}(), dc)

function delete!{K,V}(p::DictChain{K,V}, key::K)
  for d in p.dicts
    delete!(d, key)
  end
  p
end

# The default parameters storage in BlackBoxOptim
typealias Parameters Associative{Symbol,Any}

typealias ParamsDict Dict{Symbol,Any}
typealias ParamsDictChain DictChain{Symbol,Any}

# Default place for parameters argument in many methods
const EMPTY_PARAMS = @compat(Dict{Symbol, Any})

# FIXME Clean this up
function convert_to_dict_symbol_any(parameters)
  if isa(parameters, Dict{Any,Any}) || isa(parameters, Array{Any,1})
    # convert into Dict{Symbol, Any}
    parameters = Dict(@compat(Tuple{Symbol,Any})[(k::Symbol, v) for (k,v) in parameters])
  end
  return parameters
end
