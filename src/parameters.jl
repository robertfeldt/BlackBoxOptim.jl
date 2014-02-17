# An ordered set of dicts which we traverse to find the parameter value.
# Returns nothing if no param setting is found in any of the dicts.
type Parameters
  dicts  # First hash is searched first and then in order until the last

  Parameters(dicts...) = begin
    # Ensure we have at least one hash there if none are given
    if length(dicts) == 0
      dicts = [Dict{Any, Any}()]
    end

    new(dicts)
  end
end

function getindex(p::Parameters, key)
  stringkey = string(key)
  symbolkey = symbol(key)
  for d in p.dicts
    if haskey(d, symbolkey)
      return getindex(d, symbolkey)
    elseif haskey(d, stringkey)
      return getindex(d, stringkey)
    end
  end
  return nothing
end

function setindex!(p::Parameters, value, key)
  setindex!(first(p.dicts), value, key)
end

import Base.haskey, Base.get, Base.merge, Base.delete!

function haskey(p::Parameters, key)
  for d in p.dicts
    if haskey(d, key)
      return true
    end
  end
  return false
end

function get(p::Parameters, key, default = nothing)
  value = getindex(p, key)
  return (value == nothing ? default : value)
end

# In a merge the last parameter should be prioritized since this is the way
# the normal Julia merge of dicts works.
merge(p1::Union(Dict, Parameters), p2::Union(Dict, Parameters)) = Parameters(p2, p1)

function delete!(p::Parameters, key)
  for d in p.dicts
    delete!(d, key)
  end
end