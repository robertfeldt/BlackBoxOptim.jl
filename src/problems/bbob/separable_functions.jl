# Computational cores of BBOB separable functions f1-f5

# f1 = Sphere function
function sphere_function{T <: Number}(x::Vector{T})
    sum(x.^2, 1)
end
