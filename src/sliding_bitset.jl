"""
BitSet-based container that could be used to e.g. store the ids of completed
jobs by `MultithreadEvaluator`.
Assumes there is `max_seq_el` element, so that all *1:max_seq_el* elements
are in the set.
"""
mutable struct SlidingBitset
    max_seq_el::Int # the set contains all elements from 1 to max_seq_el
    max_el::Int     # max element contained in the set
    els::BitSet     # elements greater than max_seq_el that are in the set

    SlidingBitset(; max_seq_el::Integer=0) = new(max_seq_el, max_seq_el, BitSet())
end

# how much the BitSet offset could differ from max_seq_el/64
# might need to be increased in massively parallel environments
# (more than SLIDE_OFFSET*64 threads)
const SLIDE_OFFSET = 8

function Base.push!(set::SlidingBitset, el::Integer)
    if el > set.max_el
        set.max_el = el
    end
    if el == set.max_seq_el+1
        # update the next sequential el
        set.max_seq_el = el
        # see if max_seq_el could be further advanced using els
        while (set.max_el > set.max_seq_el) && in(set.els, set.max_seq_el+1)
            set.max_seq_el += 1
        end
        if !isempty(set.els) && Base._div64(set.max_seq_el) > set.els.offset + SLIDE_OFFSET
            # slide els bitset
            len = length(set.els.bits)
            if len > SLIDE_OFFSET
                copyto!(set.els.bits, 1, set.els.bits, SLIDE_OFFSET+1, len-SLIDE_OFFSET)
                fill!(view(set.els.bits, max(1, len-SLIDE_OFFSET+1):len), 0)
            else
                fill!(set.els.bits, 0)
            end
            set.els.offset += SLIDE_OFFSET
        end
    else
        push!(set.els, el)
    end
    return set
end

Base.in(set::SlidingBitset, el::Integer) =
    el <= set.max_seq_el || in(set.els, el)
