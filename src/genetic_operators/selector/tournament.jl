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
    candidates = sample(1:popsize(population), n_candidates; replace=false)
    Int[tournament(sel, population,
                   candidates[(floor(Int, n_candidates*(i-1)/n_tours)+1):floor(Int, n_candidates*i/n_tours)]) for i in 1:n_tours]
end

"""
    Simulate tournament among specified `candidates`.

    Returns the index of the winner.
"""
function tournament(sel::TournamentSelector, population, candidates)
    wins = zeros(length(candidates))
    for i in eachindex(candidates)
        ifitness = fitness(population, candidates[i])
        for j in (i+1):length(candidates)
            hat = sel.hat_comp(ifitness, fitness(population, candidates[j]))
            if hat == 1 # 1st wins, gets 3 pts
                wins[i] += 3
            elseif hat == -1 # 2nd wins, gets 3 pts
                wins[j] += 3
            else # draw, each gets 1 pt
                wins[i] += 1
                wins[j] += 1
            end
        end
    end
    candidates[sortperm(wins, rev=true)[1]]
end
