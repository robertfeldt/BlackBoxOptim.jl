
#####################################################################
# Base functions.
#####################################################################
sphere(x) = sum(abs2, x)

"""
Schwefel's ellipsoid.
"""
function ellipsoid(x)
    res = 0.0
    cumsum = 0.0
    for xx in x
        cumsum += xx
        res += abs2(cumsum)
    end
    return res
end

function elliptic(x)
    D = length(x)
    condition = 1e+6
    coefficients = condition .^ range(0, stop=1, length=D)
    return sum(coefficients .* abs2.(x))
end

#= FIXME this should be more efficient implementation, but switching to it would reset the benchmarks
function elliptic(x)
    condition = 1e+6
    powers = range(0, stop=1, length=length(x))
    @inbounds return sum(i -> (condition^powers[i])*abs2(x[i]), eachindex(x))
end
=#

function rastrigin(x)
    D = length(x)
    10 * D + sum(abs2, x) - 10 * sum(xx -> cos(2π * xx), x)
end

function ackley(x)
    D = length(x)
    try
        return 20 - 20*exp(-0.2.*sqrt(sum(abs2, x)/D)) - exp(sum(xx -> cos(2π * xx), x)/D) + e
    catch ex
        # Sometimes we have gotten a DomainError from the cos function so we protect this call
        println(ex)
        println("For input x = ", x)
        # Return a very large fitness value to indicate that this is NOT the solution we want.
        # TODO: Fix better way to handle this!
        return 9.99e100
    end
end

function schwefel1_2(x)
    D = length(x)
    partsums = zeros(D)
    partsum = 0
    for i in 1:D
        partsum += x[i]
        partsums[i] = partsum
    end
    return sum(abs2, partsums)
end

function rosenbrock(x)
    n = length(x)
    return sum(100*(view(x, 2:n) .- view(x, 1:(n-1)).^2).^2 .+ (view(x, 1:(n-1)) .- 1).^2)
end

#= FIXME this should be more efficient implementation, but switching to it would reset the benchmarks
function rosenbrock(x)
    @inbounds return sum(i -> 100.0*abs2(x[i+1] - x[i]^2) + abs2(x[i] - 1.0),
                         1:length(x)-1)
end
=#

step(x) = sum(xx -> ceil(xx + 0.5), x)

function griewank(x)
    n = length(x)
    1 + (1/4000)*sum(abs2, x) - prod(cos.(x ./ sqrt.(1:n)))
end

schwefel2_22(x) = sum(abs, x) + prod(abs, x)

schwefel2_21(x) = maximum(abs, x)

# I'm unsure about this one since it does not return the expected minima at
# [1.0, 1.0].
function schwefel2_26(x)
    D = length(x)
    418.98288727243369 * D - sum(xx -> xx * sin(sqrt(abs(xx))), x)
end

cigar(x) = x[1]^2 + 1e6 * sum(abs2, view(x, 2:length(x)))

cigtab(x) = x[1]^2 + 1e8 * x[end]^2 + 1e4 * sum(abs2, view(x, 2:length(x)-1))

"""
Generic function to define `ShekelN` problems.
"""
function shekel(x, a, c)
    @assert length(x) == size(a, 2)
    @assert length(c) == size(a, 1)
    sum = 0.0
    @inbounds for i in eachindex(c)
        den = 0.0
        for j in 1:size(a, 2)
            den += abs2(x[j] - a[i,j])
        end
        sum = sum - 1 / (den + c[i])
    end
    return sum
end

# constants for `Shekel10`
Shekel10_A = [4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7; 2 9 2 9; 5 5 3 3; 8 1 8 1; 6 2 6 2; 7 3.6 7 3.6]
Shekel10_C = [0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5]

"""
`Shekel10` is a 4D, multi-minima, non-separable test problem. Our implementation
is based on the C code in:
http://www.math.ntu.edu.tw/~wwang/cola_lab/test_problems/multiple_opt/multiopt_prob/Shekel10/Shekel10.c
"""
shekel10(x) = shekel(x, Shekel10_A, Shekel10_C)

# constants for `Shekel7`
Shekel7_A = [4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7; 2 9 2 9; 5 5 3 3]
Shekel7_C = [0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3]

"""
`Shekel7` is a 4D, multi-minima, non-separable test problem. Our implementation
is based on the C code in:
http://www.math.ntu.edu.tw/~wwang/cola_lab/test_problems/multiple_opt/multiopt_prob/Shekel7/Shekel7.c
"""
shekel7(x) = shekel(x, Shekel7_A, Shekel7_C)

# constants for `Shekel5`
Shekel5_A = [4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7]
Shekel5_C = [0.1, 0.2, 0.2, 0.4, 0.4]

"""
`Shekel5` is a 4D, multi-minima, non-separable test problem. Our implementation
is based on the C code in:
http://www.math.ntu.edu.tw/~wwang/cola_lab/test_problems/multiple_opt/multiopt_prob/Shekel5/Shekel5.c
"""
shekel5(x) = shekel(x, Shekel5_A, Shekel5_C)

"""
Generic function for `Hartman N` problem family.
"""
function hartman(x, alpha, A, P)
    @assert size(A) == size(P)
    @assert length(x) == size(P, 2)
    @assert length(alpha) == size(P, 1)
    res = 0.0
    @inbounds for i in 1:length(alpha)
        isum = 0.0
        for j in 1:size(A, 2)
            isum += A[i,j] * (x[j] - P[i,j])^2
        end
        res -= alpha[i] * exp(-isum)
    end
    res
end

# constants for `Hartman6`
Hartman6_alpha = [1.0 1.2 3.0 3.2]
Hartman6_A = [10 3 17 3.50 1.7 8; 0.05 10 17 0.1 8 14; 3 3.5 1.7 10 17 8; 17 8 0.05 10 0.1 14]
Hartman6_P = 1e-4 * [1312 1696 5569 124 8283 5886; 2329 4135 8307 3736 1004 9991; 2348 1451 3522 2883 3047 6650; 4047 8828 8732 5743 1091 381]

"""
`Hartman 6D` is a multi-minima, non-separable test problem. Our implementation is based on:
http://www.sfu.ca/~ssurjano/hart6.html
"""
hartman6(x) = hartman(x, Hartman6_alpha, Hartman6_A, Hartman6_P)

# constants for `Hartman3`
Hartman3_alpha = [1.0 1.2 3.0 3.2]
Hartman3_A = [3.0 10 30; 0.1 10 35; 3.0 10 30; 0.1 10 36]
Hartman3_P = 1e-4 * [3689 1170 2673; 4699 4387 7470; 1091 8732 5547; 381 5743 8828]

"""
`Hartman 3D` is a multi-minima, non-separable test problem. Our implementation is based on:
http://www.sfu.ca/~ssurjano/hart3.html
However, we get a different global minima than the one stated on that page.
The minima should be -3.86278 but we get a different one so use that in the problem spec:
"""
hartman3(x) = hartman(x, Hartman3_alpha, Hartman3_A, Hartman3_P)

#####################################################################
# S2 functions in addition to the base functions above. As stated
# in table II of the JADE paper: http://150.214.190.154/EAMHCO/pdf/JADE.pdf
#####################################################################

quartic(x) = sum((i,x) -> i*x^4, enumerate(x))

noisy_quartic(x) = quartic(x) + rand()

s2_step(x) = sum(xx -> abs2(ceil(xx + 0.5)), x)

#####################################################################
# Misc other interesting optimization functions and families.
#####################################################################

# Schwefel 2.13 is a hard one...
#function f=schwefel_213(x)
#global initial_flag
#persistent a b A alpha
#[ps,D]=size(x);
#if initial_flag==0
#    initial_flag=1;
#    load schwefel_213_data
#    if length(alpha)>=D
#        alpha=alpha(1:D);a=a(1:D,1:D);b=b(1:D,1:D);
#    else
#        alpha=-3+6*rand(1,D);
#        a=round(-100+200.*rand(D,D));
#        b=round(-100+200.*rand(D,D));
#    end
#    alpha=repmat(alpha,D,1);
#    A=sum(a.*sin(alpha)+b.*cos(alpha),2);
#end
#for i=1:ps
#    xx=repmat(x(i,:),D,1);
#    B=sum(a.*sin(xx)+b.*cos(xx),2);
#    f(i,1)=sum((A-B).^2,1);
#end

"""
Generator for the family of deceptive functions from the
Cuccu2011 paper on novelty-based restarts. We have vectorized it to allow
more than 1D versions. The Cuccu2011 paper uses the following values for
```
(l, w) = [(5, 0),  (15, 0),  (30, 0),
          (5, 2),  (15, 2),  (30, 2),
          (5, 10), (15, 10), (30, 10)]
```
and notes that `(15, 2)` and `(30, 2)` are the most difficult instances.
"""
function deceptive_cuccu2011(l, w)
    (x) -> begin
        sumabsx = sum(abs, x)
        if sumabsx <= 1
            return sum(abs2, x)
        elseif sumabsx >= l+1
            return sum(xx -> abs2(abs(xx) - l), x)
        else
            return (1 - 0.5 * sum(xx -> abs2(sin( (π * w * (abs(xx) - 1)) / l )), x))
        end
    end
end

# Deceptive/hardest instances:
deceptive_cuccu2011_15_2 = deceptive_cuccu2011(15, 2)
deceptive_cuccu2011_30_2 = deceptive_cuccu2011(30, 2)

"""
From section 3, page 7, of Tsallis1996 paper:
    Tsallis and Stariolo, "Generalized simulated annealing", Physica A, 1996.
available from http://www.if.ufrgs.br/~stariolo/publications/TsSt96_PhysA233_395_1996.pdf
the original paper used this as a 4-dimensional problem but here it is generalized.
"""
energy_tsallis1996(x::AbstractArray) = sum(xx -> abs2(xx^2 - 8.0), x) + 5.0*sum(x) + 57.3276
