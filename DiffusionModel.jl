using JLD,Plots,DifferentialEquations,Parameters
WORKDIR = "/home/pmxjr3/NetworkDiffusionModel"
include("NetworkSetup.jl")
include("RHS.jl")
c = 13

SC,dist,lags,N,FC,missingROIs = networksetup(c;digits=3,nSC=2,nFC=1,N=140,normalise=false)
controlR2 = load("$WORKDIR/controlR2.jld","controlR2")
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
tspan = (0.,100.)
PARAMS = PARS(SC,lags,N)
prob = DDEProblem(NetworkFlow,u0,h,tspan,PARAMS)
sol = solve(prob,alg,dt=0.1)

u = sol[:,:]
