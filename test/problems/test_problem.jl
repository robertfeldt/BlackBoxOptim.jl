fsabs(x) = sum(abs(x))
fsumsq(x) = sum(x.^2)

facts("FixedDimProblem") do
  context("1-dimensional, single-objective sumabs") do
    ss = symmetric_search_space(1)
    p = BlackBoxOptim.FixedDimProblem("sumabs", [fsabs], ss, [0.0])

    @fact is_fixed_dimensional(p)         --> true
    @fact is_any_dimensional(p)           --> false
    @fact is_single_objective_problem(p)  --> true
    @fact is_multi_objective_problem(p)   --> false
    @fact numdims(p)                      --> 1
    @fact search_space(p)                 --> ss

    @fact eval1([0.0], p)                 --> 0.0
    @fact eval1([1.2], p)                 --> 1.2
    @fact eval1([-1.9], p)                --> 1.9
  end

  context("3-dimensional, single-objective sumabs") do
    ss = symmetric_search_space(3)
    p = BlackBoxOptim.FixedDimProblem("sumabs", [fsabs], ss, [0.0])

    @fact is_fixed_dimensional(p)         --> true
    @fact is_any_dimensional(p)           --> false
    @fact is_single_objective_problem(p)  --> true
    @fact is_multi_objective_problem(p)   --> false
    @fact numdims(p)                      --> 3
    @fact search_space(p)                 --> ss

    @fact eval1([0.0, 1.0, 2.0], p)       --> 3.0
    @fact eval1([-1.0, 1.0, 2.0], p)      --> 4.0
  end

  context("1-dimensional, multi-objective sumabs_sumsq") do
    ss = symmetric_search_space(1)
    p = BlackBoxOptim.FixedDimProblem("sumabs_sumsq", [fsabs, fsumsq], ss, [0.0, 0.0])

    @fact is_fixed_dimensional(p)         --> true
    @fact is_any_dimensional(p)           --> false
    @fact is_single_objective_problem(p)  --> false
    @fact is_multi_objective_problem(p)   --> true
    @fact numdims(p)                      --> 1
    @fact search_space(p)                 --> ss

    @fact eval1([0.0], p)                 --> 0.0
    @fact eval1([1.2], p)                 --> 1.2
    @fact eval1([-1.9], p)                --> 1.9

    @fact evalall([0.0], p)               --> [0.0, 0.0]
    @fact evalall([1.2], p)               --> [1.2, 1.2^2]
    @fact evalall([-1.9], p)              --> [1.9, 1.9^2]
  end

  context("Fixed dim problem from an anydim one") do
    ap = anydim_problem("sumabs", fsabs, (0.0, 1.0), 0.0)
    p1 = as_fixed_dim_problem(ap, 1)

    @fact is_fixed_dimensional(p1)         --> true
    @fact is_any_dimensional(p1)           --> false
    @fact is_single_objective_problem(p1)  --> true
    @fact is_multi_objective_problem(p1)   --> false
    @fact numdims(p1)                      --> 1

    @fact eval1([0.0], p1)                 --> 0.0
    @fact eval1([1.2], p1)                 --> 1.2
    @fact eval1([-1.9], p1)                --> 1.9

    p3 = as_fixed_dim_problem(ap, 3)

    @fact is_fixed_dimensional(p3)         --> true
    @fact is_any_dimensional(p3)           --> false
    @fact is_single_objective_problem(p3)  --> true
    @fact is_multi_objective_problem(p3)   --> false
    @fact numdims(p3)                      --> 3

    @fact eval1([0.0, 1.0, 2.0], p3)       --> 3.0
    @fact eval1([-1.0, 1.0, 2.0], p3)      --> 4.0
  end
end

facts("ShiftedAndBiasedProblem") do
  context("Only x shifted 1-dim") do
    ss = symmetric_search_space(1)
    subp = BlackBoxOptim.FixedDimProblem("sumabs", [fsabs], ss, [0.0])
    sp = BlackBoxOptim.shifted(subp)

    @fact is_fixed_dimensional(sp)         --> true
    @fact is_any_dimensional(sp)           --> false
    @fact is_single_objective_problem(sp)  --> true
    @fact is_multi_objective_problem(sp)   --> false
    @fact numdims(sp)                      --> 1

    xs = sp.xshift
    @fact eval1(xs + [0.0], sp)                 --> 0.0
    @fact eval1(xs + [1.2], sp)                 --> 1.2
    @fact eval1(xs + [-1.9], sp)                --> 1.9
  end

  context("Shifted and biased 2-dim") do
    ss = symmetric_search_space(2)
    subp = BlackBoxOptim.FixedDimProblem("sumabs", [fsabs], ss, [0.0])
    sp = BlackBoxOptim.shifted(subp; funcshift = 1.3)

    @fact is_fixed_dimensional(sp)         --> true
    @fact is_any_dimensional(sp)           --> false
    @fact is_single_objective_problem(sp)  --> true
    @fact is_multi_objective_problem(sp)   --> false
    @fact numdims(sp)                      --> 2

    xs = sp.xshift
    @fact eval1(xs + [0.0, 1.0], sp)       --> 1.0 + 1.3
    @fact eval1(xs + [1.2, -1.3], sp)      --> 2.5 + 1.3
    @fact eval1(xs + [-1.9, 0.0], sp)      --> 1.9 + 1.3
  end

  context("Within ftol") do
    ss = symmetric_search_space(2)
    subp = BlackBoxOptim.FixedDimProblem("sumabs", [fsabs], ss, [0.0])

    @fact fitness_is_within_ftol(subp, 0.1, 0.2)    --> false
    @fact fitness_is_within_ftol(subp, 0.1, 0.09)   --> true
  end
end
