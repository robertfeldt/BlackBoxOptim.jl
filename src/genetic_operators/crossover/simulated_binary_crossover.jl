"""
Simulated Binary Crossover (SBX).

See Deb&Agrawal "Simulated binary crossover for continuous search space", 1994, Complex Systems
"""
struct SimulatedBinaryCrossover <: CrossoverOperator{2,2}
    p::Float64      # probability to modify a dimension
    η::Float64      # distribution index
    η_exp::Float64  # pre-processed index

    function SimulatedBinaryCrossover(p::Number, η::Number)
        (0.0 <= p <= 1.0) || throw(ArgumentError("p must be within [0,1]"))
        η > 0.0 || throw(ArgumentError("η must be positive"))
        new(p, η, (1/(η+1)))
    end
    SimulatedBinaryCrossover(params::Parameters) =
        SimulatedBinaryCrossover(params[:SBX_p], params[:SBX_η])
end

const SBX_DefaultOptions = ParamsDict(
    :SBX_p => 0.2,
    :SBX_η => 3.0
)

# sample SPX spread factor β
function randbeta(xover::SimulatedBinaryCrossover)
    u = rand()
    return (u < 0.5 ? 2*u : 1/(2*(1-u)))^xover.η_exp
end

function apply!(xover::SimulatedBinaryCrossover,
                targets::Vector{Individual}, targetIndices::Vector{Int},
                pop, parentIndices)
    @assert length(targets) == 2
    @assert length(targetIndices) == 2
    @assert length(parentIndices) == 2
    p1ix, p2ix = parentIndices
    @inbounds for i in 1:numdims(pop)
        if rand() <= xover.p
            # crossover
            β = randbeta(xover)
            parent_avg = 0.5 * (pop[i,p1ix] + pop[i,p2ix])
            delta = 0.5 * (pop[i,p1ix] - pop[i,p2ix])
            if rand(Bool) delta = -delta end # swap the parents
            targets[1][i] = parent_avg + β * delta
            targets[2][i] = parent_avg - β * delta
        else
            # just copy the parents
            targets[1][i] = pop[i,p1ix]
            targets[2][i] = pop[i,p2ix]
        end
    end
    return targets
end
