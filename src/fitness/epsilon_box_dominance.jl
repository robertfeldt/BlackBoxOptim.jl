# Returns two booleans where the first indicates e-box dominance and the 2nd indicates
# epsilon progress
function epsilon_box_dominates_and_epsilon_progress(u, v, epsilon)
    uindex = floor(u ./ epsilon)
    vindex = floor(v ./ epsilon)
    inner_epsilon_box_dominates_and_epsilon_progress(
        u, uindex, (epsilon .* uindex),
        v, vindex, (epsilon .* vindex))
end

# FIXME
epsilon_box_dominates_and_epsilon_progress(f1::NewFitness, f2::NewFitness, epsilon) = epsilon_box_dominates_and_epsilon_progress(fitnessvalues(f1), fitnessvalues(f2), epsilon)
epsilon_box_dominates(u, v, epsilon) = epsilon_box_dominates_and_epsilon_progress(u, v, epsilon)[1]
epsilon_box_progress(u, v, epsilon) = epsilon_box_dominates_and_epsilon_progress(u, v, epsilon)[2]

# Inner boolean expression when the intermediate results have already been calculated.
# This is to speed up processing when there is a whole archive of solutions to compare
# a new solution to.
function inner_epsilon_box_dominates_and_epsilon_progress(
    u, uindex, uindextimesepsilon,
    v, vindex, vindextimesepsilon)

    if pareto_dominates_fast(uindex, vindex)
        return (true, true)
    elseif all(uindex .== vindex)
        smaller_distance = norm(u .- uindextimesepsilon) < norm(v .- vindextimesepsilon)
        return (smaller_distance, false)
    else
        return (false, false)
    end
end
