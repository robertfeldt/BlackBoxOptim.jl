mutable struct ObjectiveFuncEvaluator
  f::Function
  archive::Archive

  ObjectiveFuncEvaluator(objectiveFunc, archive = false) = begin
    archive = archive || TopListArchive()
    new(objectiveFunc, archive)
  end
end

function evaluate(candidate, evaluator::Evaluator)
  fitness = evaluator.f(candidate)
  add!(fitness, candidate, evaluator.archive)
  fitness
end