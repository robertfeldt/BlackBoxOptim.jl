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
# the normal Julia merge of dicts works.
mergeparam{K,V}(p1::Union(Dict{K,V}, DictChain{K,V}), p2::Union(Dict{K,V}, DictChain{K,V})) = DictChain(p2, p1)

function delete!{K,V}(p::DictChain{K,V}, key::K)
  for d in p.dicts
    delete!(d, key)
  end
  p
end

# The default parameters storage in BlackBoxOptim
typealias Parameters DictChain{Symbol,Any}
