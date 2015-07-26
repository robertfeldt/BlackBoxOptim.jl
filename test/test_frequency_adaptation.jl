facts("Frequency Adaptation") do

  context("returns all indices once in the first block") do
    for(reps in 1:20)
      n = rand(1:50)
      fa = BlackBoxOptim.FrequencyAdapter(n)
      block = Int64[]
      for(i in 1:n)
        mi = BlackBoxOptim.next(fa)
        @fact mi => greater_than_or_equal(1)
        @fact mi => less_than_or_equal(n)
        push!(block, mi)
        BlackBoxOptim.update!(fa, mi, rand())
      end
      @fact sort(block) => collect(1:n)
    end
  end

  context("increases the frequency of a method that has higher progress values") do
    for(reps in 1:20)
      n = rand(2:50)
      fa = BlackBoxOptim.FrequencyAdapter(n)
      up1(mi) = BlackBoxOptim.update!(fa, mi, rand() + ((mi == 1) ? 0.4 : 0.0))
      counts = zeros(Int64, n)
      uprep(rs) = begin
        for(i in 1:rs)
          mi = BlackBoxOptim.next(fa)
          counts[mi] += 1
          up1(mi)
        end
      end
      uprep(n)
      @fact counts => ones(Int64, n)
      uprep(20*n)
      @fact counts[1] => greater_than(sum(counts[2:end])/(n-1))
    end
  end
end
