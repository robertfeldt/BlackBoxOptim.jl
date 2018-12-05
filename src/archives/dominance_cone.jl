"""
Set of points in fitness space that Pareto-dominate (`MIN` is `true`) or
are dominated by (`MIN` is `false`) the given point `pt`.

The point `pt` itself is not part of its `DominanceCone`.

If the goal is to maximize each individual fitness component
(i.e. `is_minimizing(fs::FitnessScheme) == false`), the cone of
Pareto-dominating points should use `MIN = false` and Pareto-dominated points
should use `MIN = true`.
"""
struct DominanceCone{F,N,MIN} <: SI.Region{F,N}
    pt::NTuple{N,F}

    function DominanceCone{MIN}(pt::NTuple{N,F}) where {F,N,MIN}
        (MIN isa Bool) || throw(ArgumentError("MIN should be true or false (got $MIN)"))
        new{F,N,MIN}(pt)
    end
end

is_minimizing(a::DominanceCone{<:Any,<:Any,MIN}) where MIN = MIN

# function to help define contains/intersects for all MIN values
# rely on constant propagation of lt for efficient compilation
isordered(lt::Bool, a::Number, b::Number) = lt ? a < b : a > b

@generated SI.contains(a::DominanceCone{F,N,MIN}, b::NTuple{N,F}) where {F,N,MIN} =
    quote
        isordered(MIN, a.pt[1], b[1]) && return false
        anylt = isordered(MIN, b[1], a.pt[1])
        @inbounds Base.Cartesian.@nexprs $(N-1) i -> begin
            isordered(MIN, a.pt[i+1], b[i+1]) && return false
            anylt || (anylt = isordered(MIN, b[i+1], a.pt[i+1]))
        end
        return anylt
    end

SI.contains(a::DominanceCone, b::SI.Point) = SI.contains(a, b.coord)
SI.intersects(a::DominanceCone, b::SI.Point) = SI.contains(a, b)

SI.contains(a::DominanceCone, b::SI.Rect) =
    SI.contains(a, is_minimizing(a) ? b.high : b.low)
Base.in(a::SI.Region, b::DominanceCone) = SI.contains(b, a)

SI.intersects(a::DominanceCone, b::SI.Rect) =
    SI.contains(a, is_minimizing(a) ? b.low : b.high)
SI.intersects(a::SI.Rect, b::DominanceCone) = SI.intersects(b, a)
