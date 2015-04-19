using BlackBoxOptim: epsilon_dominates_clear, epsilon_dominates_fast

function test_epsilon_dominance_func(epsilon_dominates_func)
  describe("with scalar epsilon") do
    test("simple cases with float vector fitness values of length 1") do
      @check epsilon_dominates_func([0.0], [0.0], 1.0)
      @check epsilon_dominates_func([0.99], [0.0], 1.0)
  
      @check epsilon_dominates_func([-1.0], [-1.0], 1.0)
      @check epsilon_dominates_func([-1.0], [-1.0], 0.00001)
  
      @check !epsilon_dominates_func([1.0], [0.0], 1.0)
    end
  
    test("simple cases with float vector fitness values of length 2") do
      @check epsilon_dominates_func([1.0, 2.0], [1.0, 2.0], 0.10)
      @check epsilon_dominates_func([1.1, 2.0], [1.0, 2.0], 0.10)
      @check epsilon_dominates_func([0.0, 0.99], [1.0, 2.0], 1.0)
      @check epsilon_dominates_func([-1.0, 0.99], [0.0, 2.0], 1.0)
      @check epsilon_dominates_func([0.99, -10.0], [2.0, -9.0], 1.0)
  
      @check !epsilon_dominates_func([1.1, 2.1], [1.0, 2.0], 0.10)
      @check !epsilon_dominates_func([2.1, 1.1], [2.0, 1.0], 0.10)
      @check !epsilon_dominates_func([0.1, 1.1], [0.0, 1.0], 0.10)
      @check !epsilon_dominates_func([-0.9, 0.1], [-1.0, 0.0], 0.10)
      @check !epsilon_dominates_func([2.0, 3.0], [1.0, 2.0], 0.50)
    end
  
    test("simple cases with float vector fitness values of length 3") do
      @check epsilon_dominates_func([0.0, 1.0, 2.0], [1.0, 2.0, 3.0], 1.0)
  
      @check !epsilon_dominates_func([1.0, 2.0, 3.0], [0.0, 1.0, 2.0], 1.0)
  
      @check epsilon_dominates_func([0.99, 2.0,  3.0],  [0.0, 1.0, 2.0], 1.0)
      @check epsilon_dominates_func([1.0,  1.99, 3.0],  [0.0, 1.0, 2.0], 1.0)
      @check epsilon_dominates_func([1.0,  2.0,  2.99], [0.0, 1.0, 2.0], 1.0)
    end
  end

  describe("with vector epsilon") do
    test("float vector fitness values of length 2") do
      @check epsilon_dominates_func([1.0, 2.0], [1.0, 2.0], [0.1, 0.1])
      @check epsilon_dominates_func([1.09, 2.0], [1.0, 2.0], [0.1, 0.1])
      @check epsilon_dominates_func([1.1, 2.0], [1.0, 2.0], [0.1, 0.1])
      @check epsilon_dominates_func([1.1, 2.19], [1.0, 2.0], [0.1, 0.2])
      @check epsilon_dominates_func([-1.0, 2.0], [-1.0, 2.0], [0.1, 0.1])
      @check epsilon_dominates_func([-1.0, -2.0], [-1.0, -2.0], [0.1, 0.1])

      @check !epsilon_dominates_func([1.1, 2.1], [1.0, 2.0], [0.1, 0.1])
      @check !epsilon_dominates_func([1.1, 2.2], [1.0, 2.0], [0.1, 0.2])
      @check !epsilon_dominates_func([-1.0, -2.0], [-1.1, -2.5], [0.1, 0.5])
    end
  end
end

describe("epsilon_dominance") do
  describe("epsilon_dominates_clear") do
    test_epsilon_dominance_func(epsilon_dominates_clear)
  end

  describe("epsilon_dominates_fast") do
    test_epsilon_dominance_func(epsilon_dominates_fast)
  end
end