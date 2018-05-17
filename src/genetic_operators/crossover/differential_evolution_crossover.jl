abstract type DiffEvoCrossoverOperator{NP,NC} <: CrossoverOperator{NP,NC} end

# FIXME is it possible somehow to do arithmetic operations with N?
struct DiffEvoRandBin{N} <: DiffEvoCrossoverOperator{N,1}
    cr::Float64   # probability to crossover the dimension
    f::Float64    # scale parameter

    DiffEvoRandBin{N}(cr::Number, f::Number) where N = new{N}(cr, f)
    DiffEvoRandBin{N}(options::Parameters) where N = new{N}(options[:DEX_cr], options[:DEX_f])
end

const DEX_DefaultOptions = ParamsDict(
    :DEX_f => 0.6,
    :DEX_cr => 0.7
)

crossover_parameters(xover::DiffEvoRandBin, pop, target_index) = xover.cr, xover.f

const DiffEvoRandBin1 = DiffEvoRandBin{3}
const DiffEvoRandBin2 = DiffEvoRandBin{5}

function apply!(xover::DiffEvoCrossoverOperator{3,1},
                target, target_index::Int, pop, parentIndices)
    @assert length(parentIndices) == 3
    cr, f = crossover_parameters(xover, pop, target_index)
    p1ix, p2ix, p3ix = parentIndices
    # Always ensure at least one parameter is xovered
    mut_ix = rand(1:length(target))
    @inbounds for i in 1:length(target)
        if i == mut_ix || rand() <= cr
            target[i] = pop[i,p3ix] + f * (pop[i,p1ix] - pop[i,p2ix])
        elseif target_index == 0
            target[i] = pop[i,p3ix]
        end
    end
    return target
end

function apply!(xover::DiffEvoCrossoverOperator{5,1},
                target, target_index::Int, pop, parentIndices)
    @assert length(parentIndices) == 5
    cr, f = crossover_parameters(xover, pop, target_index)
    p1ix, p2ix, p3ix, p4ix, p5ix = parentIndices
    # Always ensure at least one parameter is xovered
    mut_ix = rand(1:length(target))
    @inbounds for i in 1:length(target)
        if i == mut_ix || rand() <= cr
            target[i] = pop[i,p3ix] +
                f * (pop[i,p1ix] - pop[i,p2ix]) +
                f * (pop[i,p4ix] - pop[i,p5ix])
        elseif target_index == 0
            target[i] = pop[i,p3ix]
        end
    end
    return target
end
