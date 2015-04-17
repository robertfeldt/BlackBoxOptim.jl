num_vars_to_next_mutation_point(probMutation) = int((-log(rand())) / probMutation)

# Geometric from exponential
num_vars_to_next_mutation_point2(probMutation) = iceil(log(rand()) / log(1 - probMutation))

function count_mutations_mc(pm, size = 10, numvecs = 1000, f = num_vars_to_next_mutation_point)
  nmutations = 0
  k = 1
  for i in 1:numvecs
    while k < i*size
      l = f(pm)
      nmutations += 1
      k += l
    end
  end
  nmutations
end

function count_mutations_nm(pm, size = 10, numvecs = 1000)
  nmutations = 0
  for i in 1:numvecs
    for j in 1:size
      if rand() < pm
        nmutations += 1
      end
    end
  end
  nmutations
end

using HypothesisTests

NumReps = 30
Size = 10
NumVecs = 10000

for pm in [0.10, 0.25, 0.50, 0.75, 0.90]
  mc1 = zeros(NumReps)
  mc2 = zeros(NumReps)
  nm = zeros(NumReps)
  for i in 1:NumReps
    mc1[i] = count_mutations_mc(pm, Size, NumVecs, num_vars_to_next_mutation_point)
    mc2[i] = count_mutations_mc(pm, Size, NumVecs, num_vars_to_next_mutation_point)
    nm[i] = count_mutations_nm(pm, Size, NumVecs)
  end
  pv1 = pvalue(SignTest(mc1, nm))
  pv2 = pvalue(SignTest(mc2, nm))
  @show (pm, pv1, pv2, mean(mc1 .- nm), mean(mc2 .- nm), mean(nm))
end