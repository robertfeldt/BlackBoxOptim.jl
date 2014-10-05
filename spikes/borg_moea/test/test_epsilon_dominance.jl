describe("epsilon_dominance") do
  describe("with scalar epsilon") do
    test("simple cases with float vector fitness values of length 1") do
      @check epsilon_dominates([0.0], [0.0], 1.0)
      @check epsilon_dominates([0.99], [0.0], 1.0)
  
      @check epsilon_dominates([-1.0], [-1.0], 1.0)
      @check epsilon_dominates([-1.0], [-1.0], 0.00001)
  
      @check !epsilon_dominates([1.0], [0.0], 1.0)
    end
  
    test("simple cases with float vector fitness values of length 2") do
      @check epsilon_dominates([1.0, 2.0], [1.0, 2.0], 0.10)
      @check epsilon_dominates([1.1, 2.0], [1.0, 2.0], 0.10)
      @check epsilon_dominates([0.0, 0.99], [1.0, 2.0], 1.0)
      @check epsilon_dominates([-1.0, 0.99], [0.0, 2.0], 1.0)
      @check epsilon_dominates([0.99, -10.0], [2.0, -9.0], 1.0)
  
      @check !epsilon_dominates([1.1, 2.1], [1.0, 2.0], 0.10)
      @check !epsilon_dominates([2.1, 1.1], [2.0, 1.0], 0.10)
      @check !epsilon_dominates([0.1, 1.1], [0.0, 1.0], 0.10)
      @check !epsilon_dominates([-0.9, 0.1], [-1.0, 0.0], 0.10)
      @check !epsilon_dominates([2.0, 3.0], [1.0, 2.0], 0.50)
    end
  
    test("simple cases with float vector fitness values of length 3") do
      @check epsilon_dominates([0.0, 1.0, 2.0], [1.0, 2.0, 3.0], 1.0)
  
      @check !epsilon_dominates([1.0, 2.0, 3.0], [0.0, 1.0, 2.0], 1.0)
  
      @check epsilon_dominates([0.99, 2.0,  3.0],  [0.0, 1.0, 2.0], 1.0)
      @check epsilon_dominates([1.0,  1.99, 3.0],  [0.0, 1.0, 2.0], 1.0)
      @check epsilon_dominates([1.0,  2.0,  2.99], [0.0, 1.0, 2.0], 1.0)
    end
  end

  describe("with vector epsilon") do
    test("float vector fitness values of length 2") do
      @check epsilon_dominates([1.0, 2.0], [1.0, 2.0], [0.1, 0.1])
      @check epsilon_dominates([1.09, 2.0], [1.0, 2.0], [0.1, 0.1])
      @check epsilon_dominates([1.1, 2.0], [1.0, 2.0], [0.1, 0.1])
      @check epsilon_dominates([1.1, 2.19], [1.0, 2.0], [0.1, 0.2])
      @check epsilon_dominates([-1.0, 2.0], [-1.0, 2.0], [0.1, 0.1])
      @check epsilon_dominates([-1.0, -2.0], [-1.0, -2.0], [0.1, 0.1])

      @check !epsilon_dominates([1.1, 2.1], [1.0, 2.0], [0.1, 0.1])
      @check !epsilon_dominates([1.1, 2.2], [1.0, 2.0], [0.1, 0.2])
      @check !epsilon_dominates([-1.0, -2.0], [-1.1, -2.5], [0.1, 0.5])
    end
  end
end