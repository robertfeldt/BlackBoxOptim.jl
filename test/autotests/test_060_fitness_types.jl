using BlackBoxOptim: makefitness, NewFitness, SingleObjectiveFitness, SingleNumberFitness
using BlackBoxOptim: MultiObjectiveFitness, VectorFitness
using BlackBoxOptim: fitnessvalues, characteristics, simulations, focusversion, goalversion, evaluator

describe("Fitness types") do
  @repeat test("SingleObjectiveFitness") do
    fval = randn()
    f1 = makefitness(fval)
    @check typeof(f1) <: NewFitness
    @check typeof(f1) <: SingleObjectiveFitness
    @check typeof(f1) == SingleNumberFitness{Float64}

    @check fitnessvalues(f1) == [fval]
    @check characteristics(f1) == [fval]
    @check simulations(f1) == [fval]
    @check evaluator(f1) == nothing
    @check focusversion(f1) == 0
    @check goalversion(f1) == 0
  end

  @repeat test("MultiObjectiveFitness") do
    size = rand(2:117)
    vals = randn(size)
    f2 = makefitness(vals)
    @check typeof(f2) <: NewFitness
    @check typeof(f2) <: MultiObjectiveFitness
    @check typeof(f2) == VectorFitness{Float64}    

    @check fitnessvalues(f2) == vals
    @check characteristics(f2) == vals
    @check simulations(f2) == vals
    @check evaluator(f2) == nothing
    @check focusversion(f2) == 0
    @check goalversion(f2) == 0
  end
end