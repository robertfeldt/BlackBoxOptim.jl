abstract MeshAdaptiveDirectSearch

type LTMADS <: MeshAdaptiveDirectSearch
  n::Int64
  D::Array{Int64,2}
  tau::Float64
  wminus::Int64
  wplus::Int64
  delta_m::Float64
  delta_p::Float64
  directionGenerator

  LTMADS(n, D = [eye(n,n) -eye(n,n)]) = begin
    # Default values taken from Audet2006:
    new(D, 4, -1, +1, 1, 1, LTMADSDirectionGenerator())
  end
end

# Eq 2.1 from Audet's paper:
function update_mesh_size(mads::MeshAdaptiveDirectSearch, meshSize)
  if improvedPointFound
    wk = rand(0:mads.wplus)
  else
    wk = rand(mads.wminus:-1)
  end
  meshSize * (mads.tau ^ wk)
end

# From page 203 of Audet2006:
function update_delta_m(m::LTMADS, isMinimalFrameCenter, improvedMeshPointFound)
  s = 4 # Audet2006 uses 4 but maybe we should adapt it instead?
  if isMinimalFrameCenter
    m.delta_m = m.delta_m / s     # Increases resolution
  elseif improvedMeshPointFound && (m.delta_m <= (1/s))
    m.delta_m = s * m.delta_m     # Decrease resolution
  else
    # do nothing
  end
end

# From box on page 203 in Audet2006:
type LTMADSDirectionGenerator
  cache::Dict{Int64, Array{Int64, 1}}  # Cache of previous directions used

  LTMADSDirectionGenerator() = begin
    new(Dict{Int64, Array{Int64, 1}}())
  end
end

function generate_direction(dg::LTMADSDirectionGenerator, l, n)
  if !haskey(dg.cache, l)
    # First create the general values
    twol = 2^l
    min = -twol + 1
    max = twol - 1
    bl = rand(min:max, n)

    # Then change one of the directions
    ihat = rand(1:n)
    bl[ihat] = (randbool() ? -1 : 1) * twol

    dg.cache[l] = bl
  end

  dg.cache[l]
end
