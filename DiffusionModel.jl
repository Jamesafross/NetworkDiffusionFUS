using JLD,Plots,DifferentialEquations,Parameters
include("NetworkSetup.jl")
include("RHS.jl")
c = 13000

SC,dist,lags,N,FC,missingROIs = networksetup(c;digits=3,nSC=2,nFC=1,N=140,normalise=false)


mutable struct PARS
    SC::Matrix{Float64}
    lags::Matrix{Float64}
    N::Real
end

u0 = zeros(N)
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 0.0 : u0
alg = MethodOfSteps(Euler())
tspan = (0.,100.)
PARAMS = PARS(SC,lags,N)
prob = DDEProblem(NetworkFlow,u0,h,tspan,PARAMS)
sol = solve(prob,alg,dt=0.1)

u = sol[:,:]


heatmap(u)
