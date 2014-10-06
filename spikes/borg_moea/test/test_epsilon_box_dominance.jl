describe("epsilon_box_dominance") do
  describe("with scalar epsilon") do
    test("vectors of length 2") do
      @check  epsilon_box_dominates([0.5, 0.5],  [1.5, 1.5], 1.0)
      @check  epsilon_box_dominates([0.5, 1.5],  [1.5, 1.5], 1.0)
      @check  epsilon_box_dominates([1.5, 0.5],  [1.5, 1.5], 1.0)
      @check  epsilon_box_dominates([1.1, 1.3],  [1.5, 1.5], 1.0)
      @check !epsilon_box_dominates([1.6, 1.5],  [1.5, 1.5], 1.0)
      @check !epsilon_box_dominates([1.5, 1.9],  [1.5, 1.5], 1.0)
      @check !epsilon_box_dominates([1.5, 1.5],  [1.5, 1.5], 1.0)

      @check  epsilon_box_progress([0.5, 0.5],  [1.5, 1.5], 1.0)
      @check  epsilon_box_progress([0.5, 1.5],  [1.5, 1.5], 1.0)
      @check  epsilon_box_progress([1.5, 0.5],  [1.5, 1.5], 1.0)
      @check !epsilon_box_progress([1.1, 1.3],  [1.5, 1.5], 1.0)
      @check !epsilon_box_progress([1.6, 1.5],  [1.5, 1.5], 1.0)
      @check !epsilon_box_progress([1.5, 1.9],  [1.5, 1.5], 1.0)
      @check !epsilon_box_progress([1.5, 1.5],  [1.5, 1.5], 1.0)

      @check epsilon_box_dominates([0.5, 0.5],  [1.0, 0.99], 1.0)
      @check epsilon_box_dominates([0.5, 0.75], [0.5, 1.0], 1.0)
      @check epsilon_box_dominates([0.5, 1.75], [2.5, 1.0], 1.0)

      @check !epsilon_box_dominates([0.5,  0.5],  [0.5, 0.5], 1.0)
      @check  epsilon_box_dominates([0.49, 0.5],  [0.5, 0.5], 1.0)
      @check  epsilon_box_dominates([0.5,  0.49], [0.5, 0.5], 1.0)

      @check epsilon_box_dominates([0.5, 0.5], [0.99, 0.99], 1.0)
    end

    test("vectors of length 3") do
      @check  epsilon_box_dominates([0.5, 0.5, 0.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_dominates([1.5, 0.5, 0.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_dominates([0.5, 1.5, 0.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_dominates([0.5, 0.5, 1.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_dominates([1.5, 1.5, 0.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_dominates([1.5, 0.5, 1.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_dominates([0.5, 1.5, 1.5], [1.5, 1.5, 1.5], 1.0)

      @check  epsilon_box_progress([0.5, 0.5, 0.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_progress([1.5, 0.5, 0.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_progress([0.5, 1.5, 0.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_progress([0.5, 0.5, 1.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_progress([1.5, 1.5, 0.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_progress([1.5, 0.5, 1.5], [1.5, 1.5, 1.5], 1.0)
      @check  epsilon_box_progress([0.5, 1.5, 1.5], [1.5, 1.5, 1.5], 1.0)
    end
  end

  describe("with vector epsilon") do
    test("symmetric epsilon vector") do
      @check  epsilon_box_dominates([0.5, 0.5],  [1.5, 1.5], [1.0, 1.0])
    end

    test("asymmetric epsilon vector") do
      @check  epsilon_box_dominates([0.5, 0.5],  [1.5, 1.5], [1.0, 2.0])
      @check  epsilon_box_dominates([1.5, 0.5],  [1.5, 1.5], [1.0, 2.0])

      @check  epsilon_box_progress([0.5, 0.5],  [1.5, 1.5], [1.0, 2.0])
      @check !epsilon_box_progress([1.5, 0.5],  [1.5, 1.5], [1.0, 2.0])
    end
  end
end
