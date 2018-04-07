type FitnessQueue{F}
    queue::Vector{F}
    delay::Int
    numadded::Int
    start::Int
    FitnessQueue{F}(delay::Int) = begin
        q = Array(F, delay)
        new(q, delay, 0, 1)
    end
end

function setall!{F}(q::FitnessQueue{F}, val::F)
    for i in 1:q.delay
        q.queue[i] = val
    end
end
delay{F}(q::FitnessQueue{F}) = q.delay
length{F}(q::FitnessQueue{F}) = q.delay
lastfitness{F}(q::FitnessQueue{F}) = q.queue[q.start]
nextidx{F}(q::FitnessQueue{F}) = (q.start == q.delay) ? 1 : (q.start + 1)
delayedfitness{F}(q::FitnessQueue{F}) = q.queue[nextidx(q)]
function add{F}(q::FitnessQueue{F}, val::F)
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
function greedy_or_late_acceptance_rule{F}(q::FitnessQueue{F}, f::F)
    f < delayedfitness(fq) || f < lastfitness(fq)
end

late_acceptance_rule{F}(q::FitnessQueue{F}, f::F) = f < delayedfitness(fq)

smaller_than_all{F}(q::FitnessQueue{F}, f::F) = all(v -> f < v, fq.queue)

type LateAcceptanceHillClimbing
    mutationOperators
end