abstract Archive

# A top list archive saves a top list of the best performing (best fitness)
# candidates/individuals seen.
type TopListArchive <: Archive
  size::Int64
  fitnesses::Any[]    # Top fitness values
  candidates::Any[]   # Top candidates corresponding to top fitness values

  TopListArchive(size = 10, scheme = ) = begin
    new(size, [], [])
  end
end

# Add a candidate with a fitness to the archive.
function add!(fitness, candidate, archive::Archive)
  if length(archive.fitnesses) < archive.size
    push!(fitness, archive.fitnesses)
    push!(candidates, archive.candidates)
    sort_toplist!(archive)
  elseif 

end