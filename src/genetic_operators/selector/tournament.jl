"""
Tournament selector.
"""
mutable struct TournamentSelector{H} <: IndividualsSelector
    hat_comp::H     # fitness comparison tri-valued operator
    size::Int       # tournament size
end

function TournamentSelector(fs::FitnessScheme, size::Int=2)
    h = HatCompare(fs)
    TournamentSelector{typeof(h)}(h, size)
end

# selection using `n_tours` tournaments
function select(sel::TournamentSelector, population, n_tours::Int)
    n_candidates = min(popsize(population), sel.size*n_tours)
    all_candidates = rand_indexes(1:popsize(population), n_candidates)

    res = Vector{Int}(undef, n_tours)
    tour_candidates = Vector{Int}(undef, sel.size)
    @inbounds for i in eachindex(res)
        copyto!(tour_candidates, 1, all_candidates, 1+(i-1)*sel.size, sel.size)
        res[i] = tournament(sel, population, tour_candidates)
    end
    return res
end

"""
Simulate tournament among specified `candidates`.

Returns the index of the winner.
"""
function tournament(sel::TournamentSelector, population, candidates)
    isempty(candidates) && return 0
    @inbounds begin
        winner_ix = candidates[1]
        best_fitness = fitness(population, winner_ix)
        for i in 2:length(candidates)
            cand_ix = candidates[i]
            ifitness = fitness(population, cand_ix)
            hat = sel.hat_comp(ifitness, best_fitness)
            if hat == -1
                best_fitness = ifitness
                winner_ix = cand_ix
            end
        end
        return winner_ix
    end
end
