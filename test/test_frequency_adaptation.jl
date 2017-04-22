@testset "Frequency Adaptation" begin

    @testset "returns all indices once in the first block" begin
        for reps in 1:20
            n = rand(1:50)
            fa = BlackBoxOptim.FrequencyAdapter(n)
            block = Int[]
            for i in 1:n
                mi = BlackBoxOptim.next(fa)
                @test mi >= 1
                @test mi <= n
                push!(block, mi)
                BlackBoxOptim.update!(fa, mi, rand())
            end
            @test sort(block) == collect(1:n)
        end
    end

    @testset "increases the frequency of a method that has higher progress values" begin
        for reps in 1:20
            n = rand(2:50)
            fa = BlackBoxOptim.FrequencyAdapter(n)
            up1(mi) = BlackBoxOptim.update!(fa, mi, rand() + ((mi == 1) ? 0.4 : 0.0))
            counts = zeros(Int, n)
            uprep(rs) = begin
                for i in 1:rs
                    mi = BlackBoxOptim.next(fa)
                    counts[mi] += 1
                    up1(mi)
                end
            end
            uprep(n)
            @test counts == ones(Int, n)
            uprep(20*n)
            @test counts[1] > sum(counts[2:end])/(n-1)
        end
    end
end
