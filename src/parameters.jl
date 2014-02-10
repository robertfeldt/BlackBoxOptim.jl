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
  for h in p.dicts
    if haskey(h, symbolkey)
      return getindex(h, symbolkey)
    elseif haskey(h, stringkey)
      return getindex(h, stringkey)
    end
  end
  return nothing
end

function setindex!(p::Parameters, value, key)
  setindex!(first(p.dicts), value, key)
end