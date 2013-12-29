type AdaptiveCoordinateDescent <: Optimizer
  AdaptiveCoordinateDescent(searchSpace; options = DefaultOptions) = begin
    parent = rand_individual(searchSpace)
  end
end

function setup(acd::AdaptiveCoordinateDescent, evaluator::Evaluator)
  acd.f_best = evaluate(acd.parent, evaluator)
end

function ask(acd::AdaptiveCoordinateDescent)
end

function tell(acd::AdaptiveCoordinateDescent)
end