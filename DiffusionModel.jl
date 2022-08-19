using JLD,Plots,DifferentialEquations,Parameters
WORKDIR = "/home/pmxjr3/NetworkDiffusionModel"
include("NetworkSetup.jl")
include("RHS.jl")
c = 13000

SC,dist,lags,N,FC,missingROIs = networksetup(c;digits=3,nSC=2,nFC=1,N=140,normalise=false)
controlR2 = load("$WORKDIR/controlR2.jld",controlR2)
function maxk(a, k)
    b = partialsortperm(a, 1:k, rev=true)
    return collect(zip(b, a[b]))
end

mutable struct PARS
    SC::Matrix{Float64}
    lags::Matrix{Float64}
    N::Real
end

u0 = zeros(N)
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 0.0 : u0
alg = MethodOfSteps(Euler())
tspan = (0.,2000.)
PARAMS = PARS(SC,lags,N)
prob = DDEProblem(NetworkFlow,u0,h,tspan,PARAMS)
sol = solve(prob,alg,dt=0.1)

u = sol[:,:]

sumAll = sum(u,dims=2)
sumAll[39] = 0 

k = 14

sumAllMostImpacted = maxk(sumAll[:], k)
mostImpactedModel = zeros(k)
for i = 1:k
    mostImpactedModel[i] = sumAllMostImpacted[i][1]
end

mostImpactedModel = sort(mostImpactedModel)
cyrilMostImpacted = [53, 20, 60, 54, 43, 67, 137, 50, 121, 51, 58, 21, 66, 26]
cyrilMostImpacted = sort(cyrilMostImpacted)

count = 0
for i = 1:k
    for j = 1:k
        if mostImpactedModel[j] == cyrilMostImpacted[i]
            count += 1
        end
    end
end

print("similarity = ",count,"/",k)
u[39,:] .= 0 
heatmap(u)