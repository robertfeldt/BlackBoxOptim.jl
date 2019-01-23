mutable struct FitnessQueue{F}
    queue::Vector{F}
    delay::Int
    numadded::Int
    start::Int
    FitnessQueue{F_}(delay::Int) where {F_, F} = begin
        q = Array(F, delay)
        new(q, delay, 0, 1)
    end
end

function setall!(q::FitnessQueue{F}, val::F) where F
    for i in 1:q.delay
        q.queue[i] = val
    end
end
delay(q::FitnessQueue{F}) where {F} = q.delay
length(q::FitnessQueue{F}) where {F} = q.delay
lastfitness(q::FitnessQueue{F}) where {F} = q.queue[q.start]
nextidx(q::FitnessQueue{F}) where {F} = (q.start == q.delay) ? 1 : (q.start + 1)
delayedfitness(q::FitnessQueue{F}) where {F} = q.queue[nextidx(q)]
function add(q::FitnessQueue{F}, val::F) where F
    q.numadded += 1
    if q.numadded == 1
        setall!(q, val)
        q.minval = val
    else
        nexti = nextidx(q)
        q.start = nexti
        q.queue[nexti] = val
        q.minval = min(q.minval, val)
    end
    return val
end

#fq = FitnessQueue{Float64}(3, Inf)
#lastfitness(fq)
#delayedfitness(fq)
#add(fq, 1.0)
#lastfitness(fq)
#delayedfitness(fq)
#add(fq, 2.0)
#lastfitness(fq)
#delayedfitness(fq)
#add(fq, 3.0)
#lastfitness(fq)
#delayedfitness(fq)
#add(fq, 4.0)
#lastfitness(fq)
#delayedfitness(fq)

# Burke and Yukov propose the greedy_or_late_acceptance_rule
function greedy_or_late_acceptance_rule(q::FitnessQueue{F}, f::F) where F
    f < delayedfitness(fq) || f < lastfitness(fq)
end

late_acceptance_rule(q::FitnessQueue{F}, f::F) where {F} = f < delayedfitness(fq)

smaller_than_all(q::FitnessQueue{F}, f::F) where {F} = all(v -> f < v, fq.queue)

mutable struct LateAcceptanceHillClimbing
    mutationOperators
end