"""
    haltonsequence(b, n)

Generate the first `n` numbers in Halton's sequence with base `b`.
"""
haltonsequence(b, n) = haltonnumber.(b, 1:n)

"""
    haltonnumber(base, index)

Generate the `n`-th Halton number in the sequence with base `b`.

# Note
    Implementation is based on the psudo code in:
        http://en.wikipedia.org/wiki/Halton_sequence
"""
function haltonnumber(base::Integer, index::Integer)
    res = 0
    f = 1 / base
    i = index

    while (i > 0)
        res += f * (i % base)
        i = floor(Integer, i / base)
        f = f / base
    end

    return res
end
