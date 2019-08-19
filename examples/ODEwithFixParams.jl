using Distributed
using Plots
addprocs(8)


@everywhere using DifferentialEquations 
@everywhere using DiffEqParamEstim
@everywhere using BlackBoxOptim

# Define function to fix parameters

function fixvalues(vs, fixspec, dims)
    newvs = zeros(Float64, dims)
    vsidx = fixidx = 1
    for i in 1:dims
        if (fixidx <= length(fixspec) && i == fixspec[fixidx][1])
            newvs[i] = fixspec[fixidx][2]
            fixidx += 1
        else
            newvs[i] = vs[vsidx]
            vsidx += 1
        end
    end
    newvs
end

# Define DAE System: Lorenz Attractor + Extra Variable to illustarte DAE syntax
function g(du,u,p,t)
  σ,ρ,β = p
  x,y,z,w = u
  du[1] = σ*(y-x)
  du[2] = x*(ρ-z) - y
  du[3] = x*y - β*z
  du[4] = w - (x + y + z)
end

# Define Limits
lower=Float64.([0, 0, 0])
upper=Float64.([20, 30, 10])

# Define Initial Conditions as well as time span
u0 = [1.0;0.0;0.0;1.0]
t = 0.0:0.05:1
tspan = (0.0,1.0)

# Generate Mass Matrox MM to transfrom ODE into DAE
MM=zeros(4,4)
MM[1,1]=MM[2,2]=MM[3,3]=1

# Generate Equation system
ODESystem = ODEFunction(g;mass_matrix=MM)

# Generate Function that creates ODE problem as function of parameters (optmization argument)
model_ode(p_) = ODEProblem(ODESystem, u0, tspan,p_) 

# Generate Function that solves system as function of parameters, saves every s interval and save only vars we want.
solve_model(mp_,s, i) = DifferentialEquations.solve(model_ode(mp_), Rodas5(),saveat=s, save_idxs=i)

# Generate Clean mock data / will add nisy mock data later
mock_data = Array(solve_model([10.0,28.0,8/3],0.05,[1,2,3,4]))#,4


# Create Loss function
function L(t,dat,w=nothing)  
    return L2Loss(t,dat;data_weight=w)
end

# Create objective function
loss_objective(mp_, dat)=build_loss_objective(model_ode(mp_), Rodas5(), L(t,dat))

# Simplify fitness/objective function to depend only on parameter we want to vary
function origfitness(args) 
  return loss_objective(args, mock_data)(args)
end

dims = 3

# In specific opt problem A there is one var (at position 1) with fixed value of 10.

problemfixed = [(1, 10.0)]

# Create objective function with fix value =10 at position 1
fitnessFix(vs) = origfitness(fixvalues(vs, problemfixed, dims))

# Generate options for optimization for problem with free parameters
opt_free = bbsetup(origfitness; Method=:separable_nes, SearchRange = (collect(zip(lower,upper))),
               NumDimensions = 3, MaxFuncEvals = 100000, workers=workers(), TargetFitness=100.0 )

# Generate options for optimization for problem with fix parameters (note dimensions are 2 instead of 3)
opt_fix = bbsetup(fitnessFix; Method=:separable_nes, SearchRange = (collect(zip(lower[2:3],upper[2:3]))),
               NumDimensions = 2, MaxFuncEvals = 100000, workers=workers())

# Call bboptim as usual and optimizing for fitnessFix and origfitness.
bbfix=bboptimize(opt_fix);
bbfree=bboptimize(opt_free);
bsfix=fixvalues(best_candidate(bbfix), problemfixed,dims) # then is our best solution to the general DAE problem with a fixed param.
bsfree=best_candidate(bbfree) # then is our best solution to the general DAE problem with full freedom.

# Note that within error bsfix and bs free are equal

(bsfree, origfitness(bsfree)), (bsfix, origfitness(bsfix))


# Plot results
pyplot()
# Shot Times to see Fitness to Data
Sol1=solve(ODEProblem(ODESystem, u0, (0,1.0), bsfix),Rodas5())
plot(Sol1, layout=(2,2), xlabel="Time", ylabel=[:x :y :z :w], label=["" "" "" "w=x+y+z"])
scatter!(t,[mock_data[i,:] for i in 1:4],layout=(2,2), label="")

# Long times to see attractor behavior
Sol2=solve(ODEProblem(ODESystem, u0, (0,100.0), bsfix),Rodas5())
#plot(Sol,layout=(2,2), labels=[:w :x :y "z=w+x+y"])
P1=plot(Sol2.t, [Sol2[i,:] for i in 1:4], layout=(2,2), xlabel="Time", ylabel=[:x :y :z :w], labels=["" "" "" "w=x+y+z"])
P2=plot(Sol2,vars=(1,2,3),xlabel="x", ylabel="y", zlabel="z", labels="")
l = @layout([a;b{0.5h}])
PF=plot(P1,P2, layout=l, size=(900,600))

