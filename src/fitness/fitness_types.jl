abstract type NewFitness end

abstract type SingleObjectiveFitness <: NewFitness end

abstract type MultiObjectiveFitness <: NewFitness end

# Some evaluators go through several stages of fitness evaluation:
#   1. a candidate is simulated -> *simulation(s)*
#   2. simulations have *characteristics* that are currently in *focus*
#   3. characteristics have *fitnessvalues* in relation to current *goals*
# A certain fitness value summarizes all of this information. By default
# this complex model is not really needed and we set the characteristics == fitnessvalues
# and the evaluatorversion (which determines both the focus and the goals that we active
# for this version of the evaluator) is a constant, i.e. the focus and goals do not change.
# To judge if the versions are current some fitness values saves a reference to the evaluator,
# but the default is not to have such a reference.
# For more advanced schemes you should override these methods but here are the defaults:

focusversion(f::NewFitness) = 0
goalversion(f::NewFitness) = 0
simulations(f::NewFitness) = characteristics(f)
characteristics(f::NewFitness) = fitnessvalues(f)
evaluator(f::NewFitness) = nothing

struct SingleNumberFitness{T <: Real} <: SingleObjectiveFitness
    fvalue::T
end

struct VectorFitness{T <: Real} <: MultiObjectiveFitness
    fvalues::Vector{T}
end

makefitness(fvalue::T) where {T <: Real} = SingleNumberFitness(fvalue)
makefitness(fvalues::Vector{T}) where {T <: Real} = VectorFitness(fvalues)

fitnessvalues(f::VectorFitness{T}) where {T <: Real} = f.fvalues
fitnessvalues(f::SingleNumberFitness{T}) where {T <: Real} = [f.fvalue]
