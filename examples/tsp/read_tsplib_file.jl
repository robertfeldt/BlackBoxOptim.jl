mutable struct TspLibProblem
    name::String
    numcities::Int
    weights::Matrix{Float64}
end

size(t::TspLibProblem) = t.numcities

function resetweights!(t::TspLibProblem, nc::Int)
    t.numcities = nc
    t.weights = zeros(Float64, nc, nc)
end

function cost(p::TspLibProblem, visitorder::Vector{Int})
    startcity = city = first(visitorder)
    cost = 0.0
    for i in 2:length(visitorder)
        nextcity = visitorder[i]
        cost += p.weights[city, nextcity]
        city = nextcity
    end
    cost += p.weights[city, startcity]
    return cost
end

function whenmatch(fn, re, s)
    m = match(re, s)
    if !isnothing(m)
        fn(m)
        return true
    end
    return false
end

function parse_edge_weights!(p, lines, lineindex)
    row = col = 1
    while !isnothing(match(r"^(\s*(\d+))+\s*$", lines[lineindex]))
        map(split(lines[lineindex], r"\s+")) do dstr
            if length(dstr) > 0
                p.weights[row, col] = p.weights[col, row] = parse(Int, strip(dstr))
                col += 1
                if col > row
                    row += 1
                    col = 1
                end
            end
        end
        lineindex += 1
    end
    return lineindex - 1
end

filename = "dantzig42.tsp"

function read_tsplib_file(filename)
    p = TspLibProblem("", 1, zeros(Float64, 1, 1))
    lines = readlines(filename)
    lineindex = 0
    while lineindex < length(lines)
        lineindex += 1
        l = lines[lineindex]
        whenmatch(m -> p.name = strip(m[1]), 
            r"\s*NAME\s*:\s*(.+)$", l) && continue
        whenmatch(m -> resetweights!(p, parse(Int, m[1])), 
            r"\s*DIMENSION\s*:\s*(\d+)$", l) && continue
        whenmatch(m -> lineindex = parse_edge_weights!(p, lines, lineindex+1), 
            r"\s*EDGE_WEIGHT_SECTION\s*$", l)
    end
    return p
end
