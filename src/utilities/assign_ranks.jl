"""
Assign ranks to the values but keeps the rank the same if the values are within
tolerance of each other.
"""
function assign_ranks_within_tolerance(values; by = (x) -> x, tolerance = 1e-5, rev = false)
    perm = sortperm(values, by = by, rev = rev)
    ranked = Any[]
    rank = 1
    prev = by(values[perm[1]])
    num_of_this_rank = 0

    for i in eachindex(perm)
        r = values[perm[i]]
        v = by(r)
        if abs(prev - v) > tolerance
            rank += num_of_this_rank
            num_of_this_rank = 1
        else
            num_of_this_rank += 1
        end
        push!(ranked, (rank, r, v))
        prev = v
    end

    return ranked
end
