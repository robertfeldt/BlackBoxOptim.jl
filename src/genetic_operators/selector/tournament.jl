"""
  Tournament selector.
"""
type TournamentSelector{H} <: IndividualsSelector
    hat_comp::H     # fitness comparison tri-valued operator
    size::Int       # tournament size

    function Base.call{FS}(::Type{TournamentSelector}, fs::FS, size::Int=2)
        h = HatCompare(fs)
        new{typeof(h)}(h, size)
    end
end

# selection using `n_tours` tournaments
function select(sel::TournamentSelector, population, n_tours::Int)
    n_candidates = min(popsize(population), sel.size*n_tours)
    all_candidates = sample(1:popsize(population), n_candidates; replace=false)

    res = Vector{Int}(n_tours)
    tour_candidates = Vector{Int}(sel.size)
    @inbounds for i in eachindex(res)
        copy!(tour_candidates, 1, all_candidates, 1+(i-1)*sel.size, sel.size)
        res[i] = tournament(sel, population, tour_candidates)
    end
    return res
end

"""
    Simulate tournament among specified `candidates`.

    Returns the index of the winner.
"""
function tournament(sel::TournamentSelector, population, candidates)
    if isempty(candidates)
        return 0
    end
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
