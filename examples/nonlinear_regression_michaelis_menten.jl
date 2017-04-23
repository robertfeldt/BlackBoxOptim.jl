function michaelis_menten_model(concentration, Vm, K)
    (Vm * concentration) / (K + concentration)
end

# MicMen data is taken from the README of the NLReg.jl package: https://github.com/dmbates/NLreg.jl
# 1st column is concentration and 2nd column is the Rate
MicMenData = [
    0.02 76;
    0.02 47;
    0.06 97;
    0.06 107;
    0.11 123;
    0.11 139;
    0.22 159;
    0.22 152;
    0.56 191;
    0.56 201;
    1.1  207;
    1.1  200
];

MicMenConcentration = MicMenData[:, 1];
MicMenRate = MicMenData[:, 2];

# Fitness function takes a vector of Vm and K and calculates the RSS
function mic_men_fitness(params)
    Vm, K = params
    yhat = Float64[michaelis_menten_model(c, Vm, K) for c in MicMenConcentration]
    sumabs2(MicMenRate .- yhat)
end

using BlackBoxOptim
result = bboptimize(mic_men_fitness; 
  SearchRange = (-1000.0, 1000.0), NumDimensions = 2, MaxSteps = 1e4)
Vm, K = best_candidate(result)
RSS = best_fitness(result)

println("NLReg.jl uses specific MicMen object and fit method to find:")
println("Vm = 212.684, K = 0.0641212, RSS = 1195.45")
println("\nUsing BlackBoxOptim.jl we find:")
println("Vm = $(Vm), K = $(K), RSS = $(RSS)")
